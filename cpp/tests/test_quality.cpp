#include "support.hpp"
#include "lex1dx/lex1dx.hpp"
#include "lex/lex.hpp"

#include <cmath>
#include <vector>

using lex1dx::cluster;
using lex1dx::Mode;
using lex1dx::Options;

int main() {
    test("REE vs exact is under 0.01% on small data (tuned and sig5, k in {4,8,16})", [&] {
        std::vector<double> xs;
        for (int i = 0; i < 400; ++i)
            xs.push_back(std::sin(0.31 * i) + 0.05 * (i % 17) + ((i % 5) - 2) * 0.2);
        for (int k : {4, 8, 16}) {
            const auto ex = lex::ckmeans_exact(xs, k);
            CHECK(ex.sse > 0);
            for (auto mode : {Mode::Tuned, Mode::Sig5}) {
                Options o; o.mode = mode;
                const auto r = cluster(xs, k, o);
                const double ree = (r.sse - ex.sse) / ex.sse * 100.0;
                CHECK_MSG(ree < 0.01, "k=" + std::to_string(k) + " ree=" + std::to_string(ree));
                CHECK_MSG(ree >= -1e-9, "SSE better than exact? ree=" + std::to_string(ree));
            }
        }
    });

    test("default M is max(4096, 32k)", [&] {
        CHECK(lex1dx::default_M(64) == 4096);
        CHECK(lex1dx::default_M(128) == 4096);
        CHECK(lex1dx::default_M(256) == 8192);
        CHECK(lex1dx::default_M(512) == 16384);
    });

    test("sig5 reduces Dprime on a spread-out sample", [&] {
        std::vector<double> xs;
        for (int i = 0; i < 1000; ++i) xs.push_back(std::exp(0.005 * i) * std::sin(i * 0.1 + 0.25));
        Options o; o.mode = Mode::Sig5;
        const auto r = cluster(xs, 16, o);
        CHECK(r.Dprime > 0);
        CHECK(r.Dprime <= r.D);
        CHECK(r.k() == 16);
    });

    test("coarse DP ignores a group's internal SSE consistently (sig5, M=4)", [&] {
        // Collapse to 4 groups makes {0.08, 0.63} one group with internal SSE. If the DP counts
        // that SSE only when the group shares a bin, it starts from {0,2,3,6} (SSE 1.2321...).
        // The optimum on the raw points is {0,3,5,6} with SSE 1.1746.
        const std::vector<double> xs{0.08, 0.63, 1.57, 3.27, 3.55, 4.66};
        for (auto mode : {Mode::Tuned, Mode::Sig5}) {
            Options o; o.mode = mode; o.coarse_groups = 4;
            const auto r = cluster(xs, 3, o);
            CHECK((r.bounds == std::vector<int>{0, 3, 5, 6}));
            CHECK(std::abs(r.sse - 1.1746) < 1e-12);
        }
    });

    return finish();
}
