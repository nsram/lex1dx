// Rcpp binding: one entry point over lex1dx::cluster. Input checks live in R/cluster.R.
#include <Rcpp.h>

#include "lex1dx/lex1dx.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

// [[Rcpp::export(.cluster)]]
Rcpp::List cluster_fit(Rcpp::NumericVector x, int k, std::string mode, double M, double tolerance) {
    lex1dx::Options o;
    o.tolerance = tolerance;
    if (mode == "tuned") o.mode = lex1dx::Mode::Tuned;
    else if (mode == "sig5") o.mode = lex1dx::Mode::Sig5;
    else throw std::invalid_argument("mode must be 'tuned' or 'sig5'");
    if (M < 0 || M != std::floor(M)) throw std::invalid_argument("M must be a whole number >= 0");
    o.coarse_groups = static_cast<std::size_t>(M);
    const lex1dx::Result r = lex1dx::cluster(x.begin(), static_cast<std::size_t>(x.size()), k, o);
    return Rcpp::List::create(
        Rcpp::Named("bounds") = r.bounds, Rcpp::Named("sizes") = r.sizes,
        Rcpp::Named("centers") = r.centers, Rcpp::Named("upper") = r.upper,
        Rcpp::Named("sse") = r.sse, Rcpp::Named("D") = r.D, Rcpp::Named("Dprime") = r.Dprime,
        Rcpp::Named("sig") = r.sig,
        Rcpp::Named("sweeps") = r.sweeps, Rcpp::Named("moves") = r.moves,
        Rcpp::Named("presorted") = r.presorted);
}

// [[Rcpp::export(.default_M)]]
int default_M_(int k) { return static_cast<int>(lex1dx::default_M(k)); }

// IEEE-754 bit pattern (hex of the uint64 encoding) to double, for the tests.
// [[Rcpp::export(.bits_to_double)]]
Rcpp::NumericVector bits_to_double_(Rcpp::CharacterVector s) {
    Rcpp::NumericVector out(s.size());
    for (R_xlen_t i = 0; i < s.size(); ++i) {
        const char* p = s[i];
        char* end = nullptr;
        const unsigned long long u = std::strtoull(p, &end, 16);
        while (end && (*end == ' ' || *end == '\t' || *end == '\r')) ++end;
        if (end == p || (end && *end != '\0')) {
            out[i] = NA_REAL;
            continue;
        }
        const std::uint64_t bits = static_cast<std::uint64_t>(u);
        double v = 0;
        std::memcpy(&v, &bits, sizeof(v));
        out[i] = v;
    }
    return out;
}
