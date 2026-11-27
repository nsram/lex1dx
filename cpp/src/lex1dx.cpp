#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("-ffp-contract=off")
#endif
#include "lex1dx/lex1dx.hpp"

#include "detail.hpp"
#include "siggroup.hpp"
#include "tuned_collapse.hpp"
#include "kmeans.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <memory>
#include <utility>

namespace lex1dx {
namespace {

using lex::detail::Atoms;

[[noreturn]] void fail(const std::string& msg) { throw std::invalid_argument(msg); }

std::string show(double x) {
    if (std::isnan(x)) return "NaN";
    if (std::isinf(x)) return x > 0 ? "Infinity" : "-Infinity";
    return std::to_string(x);
}

void check_finite(const double* xs, std::size_t n) {
    if (n == 0) fail("input is empty");
    for (std::size_t i = 0; i < n; ++i)
        if (!std::isfinite(xs[i]))
            fail("input value at index " + std::to_string(i) + " is " + show(xs[i]) +
                 "; only finite values are allowed");
}

void check_weights(const double* ws, std::size_t n) {
    if (ws == nullptr) fail("weights pointer is null");
    double W = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(ws[i]) || ws[i] < 0)
            fail("weight at index " + std::to_string(i) + " is " + show(ws[i]) +
                 "; weights must be finite and >= 0");
        W += ws[i];
    }
    if (!(W > 0) || !std::isfinite(W))
        fail("weights must sum to a positive finite total, got " + show(W));
}

int resolve_sig(const Options& options) {
    if (options.sig == -1) return options.mode == Mode::Sig5 ? 5 : 0;
    if (options.sig < 0 || options.sig > 9) {
        fail("sig must be an integer in 0..9, got " + std::to_string(options.sig) +
             ". Use sig=0 for full precision");
    }
    return options.sig;
}

int count_distinct(const std::vector<double>& sorted) {
    if (sorted.empty()) return 0;
    int d = 1;
    for (std::size_t i = 1; i < sorted.size(); ++i)
        if (sorted[i] != sorted[i - 1]) ++d;
    return d;
}

void apply_user_weights(siggroup::Groups& g, const std::vector<double>& sorted,
                        const std::vector<double>& wts) {
    const int G = static_cast<int>(g.start.size()) - 1;
    for (int i = 0; i < G; ++i) {
        const int lo = g.start[i], hi = g.start[i + 1];
        double W = 0, S = 0;
        for (int e = lo; e < hi; ++e) {
            W += wts[e];
            S += wts[e] * sorted[e];
        }
        g.w[i] = W;
        g.mean[i] = W > 0 ? S / W : g.mean[i];
    }
}

struct Agg { std::vector<double> w, wx, wx2; };

Agg agg_atoms(const Atoms& a, const std::vector<int>& groups) {
    Agg g;
    const int G = static_cast<int>(groups.size()) - 1;
    g.w.resize(G); g.wx.resize(G); g.wx2.resize(G);
    for (int i = 0; i < G; ++i) {
        const int lo = groups[i], hi = groups[i + 1];
        // Each group enters the DP as a point mass at its mean: wx2 = W * mu^2, not the group's
        // raw sum of squares. Wilber's interval cost is 0 for a single item, so keeping a group's
        // internal SSE in wx2 would count it only when the group shares a bin, which can change
        // the argmin. Dropping it everywhere subtracts the same constant from every partition.
        if (a.weighted) {
            const double W = a.wsum(lo, hi);
            const double mu = a.mean(lo, hi);
            g.w[i] = W; g.wx[i] = W * mu; g.wx2[i] = W * mu * mu;
        } else {
            const double c = hi - lo;
            const double mu = a.mean(lo, hi);
            g.w[i] = c; g.wx[i] = c * mu; g.wx2[i] = c * mu * mu;
        }
    }
    return g;
}

std::vector<int> wilber_to_bounds(const std::vector<int>& groups,
                                  const std::vector<std::size_t>& path) {
    std::vector<std::size_t> p = path;
    std::sort(p.begin(), p.end());
    p.erase(std::unique(p.begin(), p.end()), p.end());
    std::vector<int> b;
    const int G = static_cast<int>(groups.size()) - 1;
    for (std::size_t gi : p)
        if (gi <= static_cast<std::size_t>(G)) b.push_back(groups[static_cast<int>(gi)]);
    if (b.empty() || b.front() != 0) b.insert(b.begin(), 0);
    if (b.back() != groups[G]) b.push_back(groups[G]);
    return b;
}

std::vector<int> ckmeans_on_groups(const Atoms& a, const std::vector<int>& groups, int kk) {
    if (a.weighted) return lex::detail::exact_weighted(a, groups, kk);
    const int G = static_cast<int>(groups.size()) - 1;
    if (G == a.d && groups == a.start) return lex::detail::exact(a, kk);
    std::vector<double> means(static_cast<std::size_t>(G));
    for (int g = 0; g < G; ++g) means[g] = a.mean(groups[g], groups[g + 1]);
    return lex::detail::exact(groups, means, kk);
}

bool valid_k_bounds(const std::vector<int>& b, int kk) {
    if (static_cast<int>(b.size()) != kk + 1) return false;
    for (int i = 0; i < kk; ++i)
        if (!(b[static_cast<std::size_t>(i)] < b[static_cast<std::size_t>(i) + 1]))
            return false;
    return true;
}

std::vector<int> ckmeans_on_groups_guarded(const Atoms& a, const std::vector<int>& groups, int kk) {
    const int G = static_cast<int>(groups.size()) - 1;
    check_ckmeans_fallback_budget(kk, G);
    return ckmeans_on_groups(a, groups, kk);
}

// Try Wilber; on empty-interval / iteration-cap / other runtime_error, a null
// result, or a bound vector that is not a valid k-partition, fall back to
// Ckmeans. The memory cap is checked in the fallback and is not swallowed here.
// *used_fallback is set when Ckmeans actually ran.
std::vector<int> wilber_on_groups(const Atoms& a, const std::vector<int>& groups, int kk,
                                  bool* used_fallback) {
    if (used_fallback) *used_fallback = false;
    auto fallback = [&]() {
        if (used_fallback) *used_fallback = true;
        return ckmeans_on_groups_guarded(a, groups, kk);
    };
    try {
        const Agg agg = agg_atoms(a, groups);
        kmeans_wilber fn(agg.w, agg.wx, agg.wx2);
        fn.set_search_strategy(search_strategy::INTERPOLATION);
        auto res = fn.compute_and_report(static_cast<std::size_t>(kk));
        if (!res) return fallback();
        auto b = wilber_to_bounds(groups, res->path);
        if (!valid_k_bounds(b, kk)) return fallback();
        return b;
    } catch (const std::runtime_error&) {
        return fallback();
    }
}

std::vector<int> map_raw(const std::vector<int>& b, const std::vector<int>& elem_start, bool sig_grouped) {
    if (!sig_grouped) return b;
    std::vector<int> o(b.size());
    for (std::size_t c = 0; c < b.size(); ++c) o[c] = elem_start[b[c]];
    return o;
}

Result summarize_raw(const std::vector<double>& sorted, std::vector<int> b,
                     int D, int Dprime, int sig, lex::detail::RefineStats stats, bool presorted,
                     const std::string& solver_used) {
    const int k = static_cast<int>(b.size()) - 1;
    Result r;
    r.sizes.resize(k);
    r.centers.resize(k);
    r.upper.resize(k);
    // two-pass, centred, same form as lex1d
    double center = 0;
    const int n = static_cast<int>(sorted.size());
    for (double x : sorted) center += x;
    center /= n;
    double sse = 0;
    for (int i = 0; i < k; ++i) {
        const int lo = b[i], hi = b[i + 1];
        double s = 0;
        for (int e = lo; e < hi; ++e) s += sorted[e] - center;
        const double m = s / (hi - lo);
        for (int e = lo; e < hi; ++e) {
            const double dv = sorted[e] - center - m;
            sse += dv * dv;
        }
        r.sizes[i] = hi - lo;
        r.centers[i] = center + m;
        r.upper[i] = sorted[hi - 1];
    }
    r.bounds = std::move(b);
    r.sse = sse;
    r.D = D;
    r.Dprime = Dprime;
    r.sig = sig;
    r.sweeps = stats.sweeps;
    r.moves = stats.moves;
    r.presorted = presorted;
    r.solver_used = solver_used;
    return r;
}

Result summarize_weighted(const std::vector<double>& sorted, const std::vector<double>& w,
                          std::vector<int> b, int D, int Dprime, int sig,
                          lex::detail::RefineStats stats, bool presorted,
                          const std::string& solver_used) {
    const int k = static_cast<int>(b.size()) - 1;
    const int n = static_cast<int>(sorted.size());
    Result r;
    r.sizes.resize(k);
    r.centers.resize(k);
    r.upper.resize(k);
    double Wtot = 0, Stot = 0;
    for (int e = 0; e < n; ++e) {
        Wtot += w[e];
        Stot += w[e] * sorted[e];
    }
    const double center = Stot / Wtot;
    double sse = 0;
    for (int i = 0; i < k; ++i) {
        const int lo = b[i], hi = b[i + 1];
        double W = 0, S = 0;
        for (int e = lo; e < hi; ++e) {
            W += w[e];
            S += w[e] * (sorted[e] - center);
        }
        const double m = W == 0 ? 0.0 : S / W;
        if (W > 0) {
            for (int e = lo; e < hi; ++e) {
                const double dv = sorted[e] - center - m;
                sse += w[e] * dv * dv;
            }
        }
        r.sizes[i] = hi - lo;
        r.centers[i] = center + m;
        r.upper[i] = sorted[hi - 1];
    }
    r.bounds = std::move(b);
    r.sse = sse;
    r.D = D;
    r.Dprime = Dprime;
    r.sig = sig;
    r.sweeps = stats.sweeps;
    r.moves = stats.moves;
    r.presorted = presorted;
    r.solver_used = solver_used;
    return r;
}

} // namespace

Result cluster(const double* xs, std::size_t n, int k, const Options& options) {
    return cluster(xs, n, k, nullptr, options);
}

Result cluster(const double* xs, std::size_t n, int k, const double* weights, const Options& options) {
    check_finite(xs, n);
    if (k < 1) fail("k must be at least 1, got " + std::to_string(k));
    if (!(options.tolerance > 0 && options.tolerance < 1))
        fail("tolerance must be in (0, 1), got " + show(options.tolerance));
    int used_sig = resolve_sig(options);
    const bool weighted = (weights != nullptr);
    if (weighted) check_weights(weights, n);

    std::vector<double> sorted(xs, xs + n);
    std::vector<double> wts;
    bool did_sort;
    if (weighted) {
        wts.assign(weights, weights + n);
        did_sort = lex::detail::sort_pairs(sorted, wts);
    } else {
        did_sort = lex::detail::sort_values(sorted);
    }
    const int D = count_distinct(sorted);

    std::vector<int> elem_start;
    std::unique_ptr<Atoms> storage;
    bool map_groups = false;
    if (used_sig > 0) {
        // Raise precision until there are at least k groups, or we reach raw atoms.
        // This is the v0.1.0 k-preserving fallback, starting from the requested sig.
        int sig = used_sig;
        auto g = siggroup::group_fast(sorted, sig);
        if (weighted) apply_user_weights(g, sorted, wts);
        while (static_cast<int>(g.w.size()) < k &&
               static_cast<int>(g.w.size()) < D &&
               sig < 15) {
            ++sig;
            g = siggroup::group_fast(sorted, sig);
            if (weighted) apply_user_weights(g, sorted, wts);
        }
        if (static_cast<int>(g.w.size()) < k && static_cast<int>(g.w.size()) < D) {
            if (weighted) storage = std::make_unique<Atoms>(sorted, wts, 0.0);
            else storage = std::make_unique<Atoms>(sorted, 0.0);
            used_sig = 0;
        } else {
            used_sig = sig;
            elem_start = std::move(g.start);
            lex::detail::sort_pairs(g.mean, g.w);
            storage = std::make_unique<Atoms>(std::move(g.mean), std::move(g.w), 0.0);
            map_groups = true;
        }
    } else if (weighted) {
        storage = std::make_unique<Atoms>(sorted, wts, 0.0);
    } else {
        storage = std::make_unique<Atoms>(sorted, 0.0);
    }
    Atoms& a = *storage;
    a.presorted = !did_sort;
    const int Dprime = a.d;
    int kk = k;
    if (kk > a.d) kk = a.d;

    std::vector<int> groups;
    if (options.mode == Mode::Exact) {
        groups = a.start;
    } else {
        const std::size_t Mreq = options.coarse_groups == 0 ? default_M(k) : options.coarse_groups;
        const int m = static_cast<int>(std::max<std::size_t>(static_cast<std::size_t>(kk),
                                                             std::min<std::size_t>(Mreq, static_cast<std::size_t>(a.d))));
        if (m >= a.d) groups = a.start;
        else groups = tuned_collapse(a, m, nullptr);
    }

    const std::string requested = (options.solver == Solver::Ckmeans) ? "ckmeans" : "wilber";
    std::string solver_used = requested;
    std::vector<int> b;
    if (kk <= 1) {
        b = {0, a.n};
        solver_used = requested;
    } else if (kk >= static_cast<int>(groups.size()) - 1) {
        b = groups;
        solver_used = requested;
    } else if (options.solver == Solver::Ckmeans) {
        b = ckmeans_on_groups(a, groups, kk);
        solver_used = "ckmeans";
    } else {
        bool used_fallback = false;
        b = wilber_on_groups(a, groups, kk, &used_fallback);
        solver_used = used_fallback ? "ckmeans" : "wilber";
    }

    const auto stats = lex::detail::refine(a, b, options.tolerance);
    if (map_groups) b = map_raw(b, elem_start, true);
    if (weighted) return summarize_weighted(sorted, wts, std::move(b), D, Dprime, used_sig, stats, !did_sort, solver_used);
    return summarize_raw(sorted, std::move(b), D, Dprime, used_sig, stats, !did_sort, solver_used);
}

int Result::assign(double x) const {
    const auto it = std::lower_bound(upper.begin(), upper.end(), x);
    return it == upper.end() ? k() - 1 : static_cast<int>(it - upper.begin());
}

} // namespace lex1dx
