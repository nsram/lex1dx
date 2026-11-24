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

int count_distinct(const std::vector<double>& sorted) {
    if (sorted.empty()) return 0;
    int d = 1;
    for (std::size_t i = 1; i < sorted.size(); ++i)
        if (sorted[i] != sorted[i - 1]) ++d;
    return d;
}

struct Agg { std::vector<double> w, wx, wx2; };

Agg agg_atoms(const Atoms& a, const std::vector<int>& groups) {
    Agg g;
    const int G = static_cast<int>(groups.size()) - 1;
    g.w.resize(G); g.wx.resize(G); g.wx2.resize(G);
    for (int i = 0; i < G; ++i) {
        const int lo = groups[i], hi = groups[i + 1];
        if (a.weighted) {
            g.w[i] = a.wsum(lo, hi);
            g.wx[i] = a.wS(lo, hi);
            g.wx2[i] = a.wQ(lo, hi);
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

std::vector<int> map_raw(const std::vector<int>& b, const std::vector<int>& elem_start, bool sig5) {
    if (!sig5) return b;
    std::vector<int> o(b.size());
    for (std::size_t c = 0; c < b.size(); ++c) o[c] = elem_start[b[c]];
    return o;
}

Result summarize_raw(const std::vector<double>& sorted, std::vector<int> b,
                     int D, int Dprime, lex::detail::RefineStats stats, bool presorted) {
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
    r.sweeps = stats.sweeps;
    r.moves = stats.moves;
    r.presorted = presorted;
    return r;
}

} // namespace

Result cluster(const double* xs, std::size_t n, int k, const Options& options) {
    check_finite(xs, n);
    if (k < 1) fail("k must be at least 1, got " + std::to_string(k));
    if (!(options.tolerance > 0 && options.tolerance < 1))
        fail("tolerance must be in (0, 1), got " + show(options.tolerance));

    std::vector<double> sorted(xs, xs + n);
    const bool did_sort = lex::detail::sort_values(sorted);
    const int D = count_distinct(sorted);
    const bool sig5 = (options.mode == Mode::Sig5);

    std::vector<int> elem_start;
    std::unique_ptr<Atoms> storage;
    if (sig5) {
        auto g = siggroup::group_fast(sorted, 5);
        elem_start = std::move(g.start);
        lex::detail::sort_pairs(g.mean, g.w);
        storage = std::make_unique<Atoms>(std::move(g.mean), std::move(g.w), 0.0);
    } else {
        storage = std::make_unique<Atoms>(sorted, 0.0);
    }
    Atoms& a = *storage;
    a.presorted = !did_sort;
    const int Dprime = a.d;
    int kk = k;
    if (kk > a.d) kk = a.d;

    const std::size_t Mreq = options.coarse_groups == 0 ? default_M(k) : options.coarse_groups;
    const int m = static_cast<int>(std::max<std::size_t>(static_cast<std::size_t>(kk),
                                                         std::min<std::size_t>(Mreq, static_cast<std::size_t>(a.d))));

    std::vector<int> groups;
    if (m >= a.d) groups = a.start;
    else groups = tuned_collapse(a, m, nullptr);

    std::vector<int> b;
    if (kk <= 1) {
        b = {0, a.n};
    } else if (kk >= static_cast<int>(groups.size()) - 1) {
        b = groups;
    } else {
        const Agg agg = agg_atoms(a, groups);
        kmeans_wilber fn(agg.w, agg.wx, agg.wx2);
        fn.set_search_strategy(search_strategy::INTERPOLATION);
        auto res = fn.compute_and_report(static_cast<std::size_t>(kk));
        if (!res) fail("weighted wilber DP failed");
        b = wilber_to_bounds(groups, res->path);
    }

    const auto stats = lex::detail::refine(a, b, options.tolerance);
    return summarize_raw(sorted, map_raw(b, elem_start, sig5), D, Dprime, stats, !did_sort);
}

int Result::assign(double x) const {
    const auto it = std::lower_bound(upper.begin(), upper.end(), x);
    return it == upper.end() ? k() - 1 : static_cast<int>(it - upper.begin());
}

} // namespace lex1dx
