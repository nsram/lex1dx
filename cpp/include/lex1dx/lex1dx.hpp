// lex1dx: fast 1-D k-means with a 0.01% relative-excess-SSE tolerance.
// Defaults differ from lex1d (the paper-faithful exact-quality reference).
#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace lex1dx {

/// Tuned: full-precision atoms unless Options::sig is set. Sig5: alias for sig = 5.
/// Exact: no Collapse. The chosen Solver runs on the atoms (after optional grouping).
enum class Mode { Tuned, Sig5, Exact };

/// Wilber: modified Grønlund wilber-interp DP (the default). Ckmeans: Fisher / ckmeans_exact DP.
enum class Solver { Wilber, Ckmeans };

struct Options {
    Mode mode = Mode::Tuned;
    /// Significant digits for grouping. -1 means follow mode: 5 if mode is Sig5, otherwise 0.
    /// 0 is full precision (no grouping). 1..9 rounds to that many significant digits.
    /// If that grouping yields fewer than k groups and the raw data has at least k
    /// distinct values, precision is raised (up to 15, then raw atoms) so k is kept.
    int sig = -1;
    /// A move is applied only if it lowers SSE by more than tolerance times the size of its terms.
    double tolerance = 1e-11;
    /// Collapse target M. Zero means the default max(4096, 32 * k). Ignored when mode is Exact.
    std::size_t coarse_groups = 0;
    Solver solver = Solver::Wilber;
};

/// Default Collapse target: M = max(4096, 32 * k).
inline std::size_t default_M(int k) {
    if (k < 1) return 4096;
    const auto kk = static_cast<std::size_t>(k);
    const auto thirty_two_k = kk > (static_cast<std::size_t>(-1) / 32) ? static_cast<std::size_t>(-1) : 32 * kk;
    return thirty_two_k < 4096 ? 4096 : thirty_two_k;
}

/// A fitted partition of the sorted data. Bin i holds sorted values bounds[i] .. bounds[i+1]-1.
struct Result {
    std::vector<int> bounds;      ///< k + 1 cut indices into the sorted data; bounds[0] = 0, bounds[k] = n
    std::vector<int> sizes;       ///< values per bin (element counts, not sums of weights)
    std::vector<double> centers;  ///< bin means of the raw values (weighted means if weights were given)
    std::vector<double> upper;    ///< largest raw value in each bin
    double sse = 0;               ///< sum of squared deviations from the bin means (weighted if weights were given)
    int D = 0;                    ///< distinct values after sort and -0 canonicalization
    int Dprime = 0;               ///< atoms after the grouping actually used (equals D when sig is 0)
    int sig = 0;                  ///< significant digits actually used; 0 means full precision
                                  ///< (requested sig was 0, or the fallback reached raw atoms)
    int sweeps = 0;               ///< refinement sweeps
    int moves = 0;                ///< refinement moves applied
    bool presorted = false;       ///< the input was already in ascending order
    std::string solver_used;      ///< "wilber" or "ckmeans": the DP that ran, or the requested
                                  ///< solver if the DP was skipped (k <= 1 or k >= number of groups)

    int k() const { return static_cast<int>(sizes.size()); } ///< bins used: min(requested k, atoms)
    /// Bin index for x: the first bin whose upper value is >= x (the last bin if none).
    int assign(double x) const;
};

/// Cluster `xs` into `k` contiguous bins of the sorted values.
/// Pipeline: sort, optional significant-digit grouping, tuned batch Collapse
/// (0.75-quantile batches, heap finish at 16384 remaining merges), a DP on the
/// M groups (wilber-interp by default, or the Fisher / ckmeans DP), then Refine.
/// mode Exact skips Collapse and runs the chosen solver on the atoms.
/// k larger than the number of atoms is clamped to it.
/// Throws std::invalid_argument on empty input, non-finite values, or k < 1.
/// A Wilber-to-Ckmeans fallback that would allocate a back-pointer table over
/// 2 GiB throws std::runtime_error instead of allocating.
Result cluster(const double* xs, std::size_t n, int k, const Options& options = {});
inline Result cluster(const std::vector<double>& xs, int k, const Options& options = {}) {
    return cluster(xs.data(), xs.size(), k, options);
}

/// Weighted variants. `weights` is a non-negative finite weight per value (same length as `xs`).
/// A null weight pointer, or the overloads above, keeps the unweighted solver.
/// `Result::sse` is then the weighted SSE, and `Result::centers` are the weighted bin means.
/// `Result::sizes` remain element counts. Weights that are negative, non-finite, all zero, or
/// of the wrong length are rejected.
Result cluster(const double* xs, std::size_t n, int k, const double* weights, const Options& options);
inline Result cluster(const std::vector<double>& xs, const std::vector<double>& weights, int k,
                      const Options& options = {}) {
    if (weights.size() != xs.size()) throw std::invalid_argument("weights size must match input");
    return cluster(xs.data(), xs.size(), k, weights.data(), options);
}

/// Bytes for the Ckmeans (Fisher DP) back-pointer table: (k-1)*(G+1)*sizeof(int).
/// Overflow saturates to UINT64_MAX. k <= 1 or G < 0 yields 0.
inline std::uint64_t ckmeans_table_bytes(int k, int G) {
    if (k <= 1 || G < 0) return 0;
    const std::uint64_t rows = static_cast<std::uint64_t>(k - 1);
    const std::uint64_t cols = static_cast<std::uint64_t>(G) + 1;
    if (cols != 0 && rows > (std::numeric_limits<std::uint64_t>::max() / cols))
        return std::numeric_limits<std::uint64_t>::max();
    const std::uint64_t cells = rows * cols;
    const std::uint64_t sz = static_cast<std::uint64_t>(sizeof(int));
    if (sz != 0 && cells > (std::numeric_limits<std::uint64_t>::max() / sz))
        return std::numeric_limits<std::uint64_t>::max();
    return cells * sz;
}

/// Scala unweighted Exact path: (k+1)*(G+1)*sizeof(int) (two unused rows).
inline std::uint64_t ckmeans_table_bytes_scala_unweighted(int k, int G) {
    if (k < 0 || G < 0) return 0;
    const std::uint64_t rows = static_cast<std::uint64_t>(k) + 1;
    const std::uint64_t cols = static_cast<std::uint64_t>(G) + 1;
    if (cols != 0 && rows > (std::numeric_limits<std::uint64_t>::max() / cols))
        return std::numeric_limits<std::uint64_t>::max();
    const std::uint64_t cells = rows * cols;
    const std::uint64_t sz = static_cast<std::uint64_t>(sizeof(int));
    if (sz != 0 && cells > (std::numeric_limits<std::uint64_t>::max() / sz))
        return std::numeric_limits<std::uint64_t>::max();
    return cells * sz;
}

/// Cap on a Wilber-to-Ckmeans fallback table: 2 GiB.
inline constexpr std::uint64_t ckmeans_fallback_cap_bytes = 2ull * 1024ull * 1024ull * 1024ull;

/// Throws std::runtime_error if ckmeans_table_bytes(k, G) exceeds the 2 GiB cap.
/// Does not allocate the table.
inline void check_ckmeans_fallback_budget(int k, int G) {
    const std::uint64_t bytes = ckmeans_table_bytes(k, G);
    if (bytes <= ckmeans_fallback_cap_bytes) return;
    const double gib = static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
    throw std::runtime_error(
        std::string("Wilber interpolation search failed; the Ckmeans fallback would allocate a (k-1)x(G+1) ")
        + "back-pointer table of " + std::to_string(gib)
        + " GiB, which exceeds the 2 GiB cap. Use mode=\"tuned\", or solver=\"ckmeans\" with a smaller k."
        + " (k=" + std::to_string(k) + ", G=" + std::to_string(G) + ")");
}

} // namespace lex1dx
