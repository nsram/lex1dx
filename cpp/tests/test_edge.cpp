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

    return finish();
}
