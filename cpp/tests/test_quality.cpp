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

    return finish();
}
