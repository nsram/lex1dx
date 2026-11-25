// lex1dx: fast 1-D k-means with a 0.01% relative-excess-SSE tolerance.
// Defaults differ from lex1d (the paper-faithful exact-quality reference).
#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace lex1dx {

/// Tuned: full-precision atoms (exact ties only). Sig5: group to 5 significant digits first.
enum class Mode { Tuned, Sig5 };

struct Options {
    Mode mode = Mode::Tuned;
    /// A move is applied only if it lowers SSE by more than tolerance times the size of its terms.
    double tolerance = 1e-11;
    /// Collapse target M. Zero means the default max(4096, 32 * k).
    std::size_t coarse_groups = 0;
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
    std::vector<int> sizes;       ///< values per bin
    std::vector<double> centers;  ///< bin means of the raw values
    std::vector<double> upper;    ///< largest raw value in each bin
    double sse = 0;               ///< sum of squared deviations from the bin means (raw values)
    int D = 0;                    ///< distinct values after sort and -0 canonicalization
    int Dprime = 0;               ///< atoms after the grouping actually used (equals D for Tuned,
                                  ///< or for Sig5 after a fallback to raw distinct values)
    int sig = 0;                  ///< significant digits used for grouping; 0 means full precision
                                  ///< (Tuned, or Sig5 fallback to raw distinct-value atoms)
    int sweeps = 0;               ///< refinement sweeps
    int moves = 0;                ///< refinement moves applied
    bool presorted = false;       ///< the input was already in ascending order

    int k() const { return static_cast<int>(sizes.size()); } ///< bins used: min(requested k, atoms)
    /// Bin index for x: the first bin whose upper value is >= x (the last bin if none).
    int assign(double x) const;
};

/// Cluster `xs` into `k` contiguous bins of the sorted values.
/// Pipeline: sort, optional significant-digit grouping, tuned batch Collapse
/// (0.75-quantile batches, heap finish at 16384 remaining merges), weighted
/// wilber-interp DP on the M groups, then Refine on those atoms (sig5 then
/// maps the refined bounds back to raw indices).
/// Sig5 starts at 5 significant digits. If that yields fewer than k groups and
/// the raw data has at least k distinct values, grouping is retried at 6, 7, ...
/// digits, and if 15 digits are still not enough, raw distinct-value atoms are
/// used (the same atoms as Tuned) so the result has min(k, D) bins.
/// k larger than the number of distinct values is clamped to it.
/// Throws std::invalid_argument on empty input, non-finite values, or k < 1.
Result cluster(const double* xs, std::size_t n, int k, const Options& options = {});
inline Result cluster(const std::vector<double>& xs, int k, const Options& options = {}) {
    return cluster(xs.data(), xs.size(), k, options);
}

} // namespace lex1dx
