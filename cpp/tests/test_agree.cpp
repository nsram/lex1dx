// C++ side of the C++/Scala agreement: dumps cluster results for a fixed set
// next to the test binary as agree_cpp.txt (and per-case files).
#include "support.hpp"
#include "lex1dx/lex1dx.hpp"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

using lex1dx::cluster;
using lex1dx::Mode;
using lex1dx::Options;

static std::vector<double> load_f64(const std::string& path) {
    FILE* fp = std::fopen(path.c_str(), "rb");
    if (!fp) return {};
    std::fseek(fp, 0, SEEK_END);
    const long bytes = std::ftell(fp);
    std::fseek(fp, 0, SEEK_SET);
    if (bytes <= 0 || bytes % 8 != 0) { std::fclose(fp); return {}; }
    const std::size_t n = static_cast<std::size_t>(bytes) / 8;
    std::vector<double> x(n);
    if (std::fread(x.data(), 8, n, fp) != n) { std::fclose(fp); return {}; }
    std::fclose(fp);
    return x;
}

static std::string format_result(const std::string& name, const std::string& mode, int kreq,
                                 const lex1dx::Result& r) {
    std::ostringstream o;
    o << "case=" << name << " mode=" << mode << " k=" << kreq
      << " k_eff=" << r.k() << " D=" << r.D << " Dprime=" << r.Dprime
      << " sse=" << std::setprecision(17) << r.sse << "\nbounds=";
    for (std::size_t i = 0; i < r.bounds.size(); ++i) o << (i ? "," : "") << r.bounds[i];
    o << "\ncenters=";
    for (std::size_t i = 0; i < r.centers.size(); ++i) {
        o << (i ? "," : "");
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.17g", r.centers[i]);
        o << buf;
    }
    o << "\n";
    return o.str();
}

// A small built-in set that CI always has (no external files).
static std::vector<double> builtin() {
    std::vector<double> xs;
    for (int i = 0; i < 800; ++i)
        xs.push_back(std::sin(0.31 * i) + 0.05 * (i % 17) + ((i % 5) - 2) * 0.2);
    return xs;
}

int main() {
    const char* outdir = std::getenv("LEX1DX_AGREE_DIR");
    const std::string dir = outdir ? outdir : ".";
    std::ofstream all(dir + "/agree_cpp.txt");
    CHECK(all.good());

    auto run = [&](const std::string& name, const std::vector<double>& xs, int k, Mode mode) {
        Options o; o.mode = mode;
        const auto r = cluster(xs, k, o);
        const char* ms = mode == Mode::Sig5 ? "sig5" : "tuned";
        const std::string block = format_result(name, ms, k, r);
        all << block << "---\n";
        return r;
    };

    test("builtin tuned/sig5 dumps", [&] {
        const auto xs = builtin();
        for (int k : {4, 8, 16, 64}) {
            const auto t = run("builtin", xs, k, Mode::Tuned);
            const auto s = run("builtin", xs, k, Mode::Sig5);
            CHECK(t.k() == k);
            CHECK(s.k() == k);
        }
    });

    // Optional: synthetic 2^16 from gen_synth, and an audio slice, if the files exist.
    const char* data = std::getenv("LEX1DX_DATA");
    const std::string data_root = data ? data : "";
    if (!data_root.empty()) {
        test("optional synth 2^16 and audio slice, if present", [&] {
            int found = 0;
            for (const char* fam : {"uniform", "normal", "lognormal"}) {
                const auto xs = load_f64(data_root + "/synth/" + std::string(fam) + "_n65536.f64");
                if (xs.size() != 65536) continue;
                ++found;
                for (int k : {64, 256, 512})
                    for (auto mode : {Mode::Tuned, Mode::Sig5})
                        run(std::string(fam) + "_n65536", xs, k, mode);
            }
            auto audio = load_f64(data_root + "/audio_s0.raw.f64");
            if (audio.size() >= 65536) {
                ++found;
                audio.resize(65536);
                for (int k : {64, 256, 512})
                    for (auto mode : {Mode::Tuned, Mode::Sig5})
                        run("audio_slice65536", audio, k, mode);
            }
            CHECK_MSG(found > 0, "LEX1DX_DATA set but no files found");
        });
    }

    all.close();
    return finish();
}
