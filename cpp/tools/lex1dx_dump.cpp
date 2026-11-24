// Print bounds, centers and SSE for one clustering, in a stable text format.
#include "lex1dx/lex1dx.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

static std::vector<double> load_f64(const char* path) {
    FILE* fp = std::fopen(path, "rb");
    if (!fp) { std::fprintf(stderr, "cannot read %s\n", path); std::exit(2); }
    std::fseek(fp, 0, SEEK_END);
    const long bytes = std::ftell(fp);
    std::fseek(fp, 0, SEEK_SET);
    if (bytes <= 0 || bytes % 8 != 0) { std::fprintf(stderr, "bad f64 %s\n", path); std::exit(2); }
    const std::size_t n = static_cast<std::size_t>(bytes) / 8;
    std::vector<double> x(n);
    if (std::fread(x.data(), 8, n, fp) != n) { std::fprintf(stderr, "short read %s\n", path); std::exit(2); }
    std::fclose(fp);
    return x;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: lex1dx_dump {tuned|sig5} IN.f64 K\n");
        return 2;
    }
    const std::string mode = argv[1];
    lex1dx::Options o;
    if (mode == "sig5") o.mode = lex1dx::Mode::Sig5;
    else if (mode == "tuned") o.mode = lex1dx::Mode::Tuned;
    else { std::fprintf(stderr, "mode must be tuned or sig5\n"); return 2; }
    const int k = std::atoi(argv[3]);
    const auto xs = load_f64(argv[2]);
    const auto r = lex1dx::cluster(xs, k, o);
    std::printf("n=%zu k=%d mode=%s D=%d Dprime=%d sse=%.17g sweeps=%d moves=%d\n",
                xs.size(), r.k(), mode.c_str(), r.D, r.Dprime, r.sse, r.sweeps, r.moves);
    std::printf("bounds=");
    for (std::size_t i = 0; i < r.bounds.size(); ++i)
        std::printf("%s%d", i ? "," : "", r.bounds[i]);
    std::printf("\ncenters=");
    for (std::size_t i = 0; i < r.centers.size(); ++i)
        std::printf("%s%.17g", i ? "," : "", r.centers[i]);
    std::printf("\n");
    return 0;
}
