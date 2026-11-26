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
        }
    });

    return finish();
}
