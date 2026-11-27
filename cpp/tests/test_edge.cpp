#include "support.hpp"
#include "lex1dx/lex1dx.hpp"
#include "lex/lex.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using lex1dx::cluster;
using lex1dx::Mode;
using lex1dx::Options;
using lex1dx::Solver;

static void rejects(const std::vector<double>& xs, int k, const std::string& fragment, Options o = {}) {
    try {
        cluster(xs, k, o);
        CHECK_MSG(false, "no exception for " + fragment);
    } catch (const std::invalid_argument& e) {
        CHECK_MSG(std::string(e.what()).find(fragment) != std::string::npos, e.what());
    }
}

int main() {
    const double inf = std::numeric_limits<double>::infinity();

    test("empty input is rejected", [&] { rejects({}, 1, "empty"); });

    test("NaN and inf are rejected", [&] {
        rejects({1.0, std::nan("")}, 1, "NaN");
        rejects({inf}, 1, "Infinity");
        rejects({0.0, -inf}, 1, "Infinity");
    });

    test("k < 1 is rejected", [&] { rejects({1.0, 2.0}, 0, "k must be at least 1"); });

    test("single value", [&] {
        for (auto mode : {Mode::Tuned, Mode::Sig5}) {
            Options o; o.mode = mode;
            const auto r = cluster(std::vector<double>{3.5}, 1, o);
            CHECK(r.k() == 1);
            CHECK(r.bounds == std::vector<int>({0, 1}));
            CHECK(r.centers[0] == 3.5);
            CHECK(r.sse == 0.0);
            CHECK(r.assign(3.5) == 0);
        }
    });

    test("all ties: k is clamped to 1", [&] {
        const std::vector<double> xs(20, 2.0);
        for (auto mode : {Mode::Tuned, Mode::Sig5}) {
            Options o; o.mode = mode;
            const auto r = cluster(xs, 8, o);
            CHECK(r.k() == 1);
            CHECK(r.bounds == std::vector<int>({0, 20}));
            CHECK(r.centers[0] == 2.0);
            CHECK(r.sse == 0.0);
            CHECK(r.D == 1);
            CHECK(r.Dprime == 1);
        }
    });

    test("k >= D is clamped to D", [&] {
        const std::vector<double> xs{1.0, 2.0, 3.0, 1.0, 2.0};
        const auto r = cluster(xs, 10, {});
        CHECK(r.k() == 3);
        CHECK(r.D == 3);
        CHECK(r.sse == 0.0);
        CHECK(r.bounds.size() == 4);
    });

    test("round-trip: sorted input matches shuffled input; assign covers the range", [&] {
        std::vector<double> xs;
        for (int i = 0; i < 200; ++i) xs.push_back(std::sin(0.17 * i) * 10.0 + 0.01 * i);
        const auto r = cluster(xs, 8, {});
        std::vector<double> sorted = xs;
        std::sort(sorted.begin(), sorted.end());
        const auto q = cluster(sorted, 8, {});
        CHECK(r.bounds == q.bounds);
        CHECK(r.centers == q.centers);
        CHECK(r.sse == q.sse);
        CHECK(q.presorted);
        CHECK(!r.presorted);
        for (double x : xs) {
            const int b = r.assign(x);
            CHECK(b >= 0 && b < r.k());
        }
        CHECK(r.assign(r.upper[0]) == 0);
        CHECK(r.assign(r.upper.back()) == r.k() - 1);
    });

    test("weights match expanding integer counts; bad weights are rejected", [&] {
        const std::vector<double> vals{1.0, 2.0, 3.0};
        const std::vector<double> wts{1.0, 3.0, 1.0};
        std::vector<double> expanded{1.0, 2.0, 2.0, 2.0, 3.0};
        const auto rw = cluster(vals, wts, 2);
        const auto re = cluster(expanded, 2);
        CHECK(rw.sse == re.sse);
        CHECK(rw.centers == re.centers);
        CHECK(static_cast<int>(rw.sizes.size()) == static_cast<int>(re.sizes.size()));
        for (std::size_t i = 0; i < rw.sizes.size(); ++i) {
            double mass = 0;
            for (int e = rw.bounds[i]; e < rw.bounds[i + 1]; ++e) mass += wts[e];
            CHECK(mass == static_cast<double>(re.sizes[i]));
        }
        try {
            cluster(vals, std::vector<double>{1.0, -1.0, 1.0}, 2);
            CHECK_MSG(false, "no exception for negative weight");
        } catch (const std::invalid_argument& e) {
            CHECK_MSG(std::string(e.what()).find("weights must be finite and >= 0") != std::string::npos, e.what());
        }
        try {
            cluster(vals, std::vector<double>{1.0, std::nan(""), 1.0}, 2);
            CHECK_MSG(false, "no exception for NaN weight");
        } catch (const std::invalid_argument& e) {
            CHECK_MSG(std::string(e.what()).find("weights must be finite and >= 0") != std::string::npos, e.what());
        }
        try {
            cluster(vals, std::vector<double>{0.0, 0.0, 0.0}, 2);
            CHECK_MSG(false, "no exception for all-zero weights");
        } catch (const std::invalid_argument& e) {
            CHECK_MSG(std::string(e.what()).find("positive finite total") != std::string::npos, e.what());
        }
        try {
            (void)cluster(vals, std::vector<double>{1.0, 1.0}, 2);
            CHECK_MSG(false, "no exception for length mismatch");
        } catch (const std::invalid_argument& e) {
            CHECK_MSG(std::string(e.what()).find("weights size must match input") != std::string::npos, e.what());
        }
    });

    test("exact solvers agree on SSE; sig 0..9 is validated", [&] {
        std::vector<double> xs;
        for (int i = 0; i < 200; ++i) xs.push_back(std::sin(0.17 * i) * 10.0 + 0.01 * i);
        Options ow; ow.mode = Mode::Exact; ow.solver = Solver::Wilber;
        Options oc; oc.mode = Mode::Exact; oc.solver = Solver::Ckmeans;
        const auto rw = cluster(xs, 8, ow);
        const auto rc = cluster(xs, 8, oc);
        const double denom = std::max(rw.sse, rc.sse);
        const double rel = denom > 0 ? std::abs(rw.sse - rc.sse) / denom : 0;
        CHECK_MSG(rel < 1e-12, "solver SSE rel=" + std::to_string(rel));
        Options bad; bad.sig = 10;
        try {
            cluster(xs, 4, bad);
            CHECK_MSG(false, "no exception for sig=10");
        } catch (const std::invalid_argument& e) {
            const std::string msg(e.what());
            CHECK_MSG(msg.find("0..9") != std::string::npos, msg);
            CHECK_MSG(msg.find("sig=0") != std::string::npos, msg);
        }
        Options s0; s0.sig = 0;
        const auto r0 = cluster(xs, 8, s0);
        const auto rt = cluster(xs, 8, {});
        CHECK(r0.sig == 0);
        CHECK(r0.bounds == rt.bounds);
        CHECK(r0.sse == rt.sse);
    });

    test("athletics20 k=5: wilber and ckmeans give the same partition", [&] {
        // time-only athletics20 values (data/athletics20.csv)
        const std::vector<double> xs = {
            9.58, 9.69, 9.69, 9.72, 9.74, 9.76, 9.76, 9.76, 9.77, 9.78,
            9.79, 9.80, 9.80, 9.82, 9.83, 9.83, 9.83, 9.83, 9.84, 9.84
        };
        for (Mode mode : {Mode::Exact, Mode::Tuned}) {
            Options ow; ow.mode = mode; ow.solver = Solver::Wilber;
            Options oc; oc.mode = mode; oc.solver = Solver::Ckmeans;
            const auto rw = cluster(xs, 5, ow);
            const auto rc = cluster(xs, 5, oc);
            CHECK_MSG(rw.bounds == rc.bounds, mode == Mode::Exact ? "exact bounds" : "tuned bounds");
            const double denom = std::max(rw.sse, rc.sse);
            const double rel = denom > 0 ? std::abs(rw.sse - rc.sse) / denom : 0;
            CHECK_MSG(rel < 1e-12, (mode == Mode::Exact ? "exact rel=" : "tuned rel=") + std::to_string(rel));
            CHECK_MSG(rw.solver_used == "wilber", "athletics wilber solver_used");
            CHECK_MSG(rc.solver_used == "ckmeans", "athletics ckmeans solver_used");
        }
    });

    auto valid_partition = [](const lex1dx::Result& r, int kk) {
        CHECK(static_cast<int>(r.bounds.size()) == kk + 1);
        CHECK(r.k() == kk);
        for (int i = 0; i < kk; ++i) {
            CHECK(r.bounds[i] < r.bounds[i + 1]);
            CHECK(r.sizes[i] > 0);
        }
    };

    test("repro {0,1,2,3} k=3: wilber falls back to ckmeans, SSE 0.5", [&] {
        const std::vector<double> xs{0.0, 1.0, 2.0, 3.0};
        const auto rw = cluster(xs, 3, {});
        Options oc; oc.solver = Solver::Ckmeans;
        const auto rc = cluster(xs, 3, oc);
        CHECK(rw.sse == 0.5);
        CHECK(rc.sse == 0.5);
        CHECK(rw.bounds == rc.bounds);
        CHECK(rw.solver_used == "ckmeans");
        CHECK(rc.solver_used == "ckmeans");
        valid_partition(rw, 3);
        Options ew; ew.mode = Mode::Exact; ew.solver = Solver::Wilber;
        const auto re = cluster(xs, 3, ew);
        CHECK(re.solver_used == "ckmeans");
        CHECK(re.sse == 0.5);
    });

    test("repro 0..255 k=27: wilber falls back and matches ckmeans", [&] {
        std::vector<double> xs(256);
        for (int i = 0; i < 256; ++i) xs[i] = static_cast<double>(i);
        const auto rw = cluster(xs, 27, {});
        Options oc; oc.solver = Solver::Ckmeans;
        const auto rc = cluster(xs, 27, oc);
        CHECK(rw.solver_used == "ckmeans");
        CHECK(rc.solver_used == "ckmeans");
        CHECK(rw.bounds == rc.bounds);
        CHECK(rw.sse == rc.sse);
        valid_partition(rw, 27);
    });

    test("skipped DP: solver_used mirrors the requested solver", [&] {
        // all ties: k is clamped to 1, DP is skipped
        const std::vector<double> ties(20, 2.0);
        const auto rw = cluster(ties, 8, {});
        Options oc; oc.solver = Solver::Ckmeans;
        const auto rc = cluster(ties, 8, oc);
        CHECK(rw.k() == 1);
        CHECK(rw.solver_used == "wilber");
        CHECK(rc.solver_used == "ckmeans");
        // k >= D: DP skipped, groups used as bins
        const std::vector<double> xs{1.0, 2.0, 3.0, 1.0, 2.0};
        const auto cw = cluster(xs, 10, {});
        Options cc; cc.solver = Solver::Ckmeans;
        const auto cck = cluster(xs, 10, cc);
        CHECK(cw.k() == 3);
        CHECK(cw.solver_used == "wilber");
        CHECK(cck.solver_used == "ckmeans");
    });

    test("ckmeans_table_bytes and fallback cap fire without allocating", [&] {
        CHECK(lex1dx::ckmeans_table_bytes(1, 100) == 0);
        CHECK(lex1dx::ckmeans_table_bytes(512, 16384) == 511ull * 16385ull * sizeof(int));
        CHECK(lex1dx::ckmeans_table_bytes(512, 16384) < lex1dx::ckmeans_fallback_cap_bytes);
        const auto big = lex1dx::ckmeans_table_bytes(512, 7000000);
        CHECK(big == 511ull * 7000001ull * sizeof(int));
        CHECK(big > lex1dx::ckmeans_fallback_cap_bytes);
        // Guard: estimate only, no table allocation. Use a tiny cap check via the public formula.
        // The runtime guard is exercised by calling cluster in Exact mode with huge D would
        // be too heavy; we test the predicate the implementation uses.
        // Guard fires without allocating the table.
        const int G_over = static_cast<int>(lex1dx::ckmeans_fallback_cap_bytes / 4); // bytes = (G+1)*4 = cap+4
        const auto over = lex1dx::ckmeans_table_bytes(2, G_over);
        CHECK(over > lex1dx::ckmeans_fallback_cap_bytes);
        lex1dx::check_ckmeans_fallback_budget(512, 16384); // under cap: no throw
        try {
            lex1dx::check_ckmeans_fallback_budget(2, G_over);
            CHECK_MSG(false, "cap guard did not throw");
        } catch (const std::runtime_error& e) {
            const std::string msg(e.what());
            CHECK_MSG(msg.find("2 GiB cap") != std::string::npos, msg);
            CHECK_MSG(msg.find("mode=\"tuned\"") != std::string::npos, msg);
        }
    });

    return finish();
}
