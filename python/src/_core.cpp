// nanobind bindings: numpy float64 arrays in, a dict of numpy arrays and scalars out.
#include "lex1dx/lex1dx.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;
using Vec = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

namespace {

template <typename T>
nb::ndarray<nb::numpy, T, nb::ndim<1>> to_numpy(const std::vector<T>& v) {
    T* data = new T[v.size() ? v.size() : 1];
    if (!v.empty()) std::memcpy(data, v.data(), v.size() * sizeof(T));
    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<T*>(p); });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(data, {v.size()}, owner);
}

nb::dict to_dict(const lex1dx::Result& r) {
    nb::dict d;
    std::vector<long long> bounds(r.bounds.begin(), r.bounds.end());
    std::vector<long long> sizes(r.sizes.begin(), r.sizes.end());
    d["bounds"] = to_numpy(bounds);
    d["sizes"] = to_numpy(sizes);
    d["centers"] = to_numpy(r.centers);
    d["upper"] = to_numpy(r.upper);
    d["sse"] = r.sse;
    d["D"] = r.D;
    d["Dprime"] = r.Dprime;
    d["sweeps"] = r.sweeps;
    d["moves"] = r.moves;
    d["presorted"] = r.presorted;
    return d;
}

nb::dict cluster_fit(Vec x, int k, const std::string& mode, std::size_t M, double tolerance) {
    lex1dx::Options o;
    o.tolerance = tolerance;
    if (mode == "tuned") o.mode = lex1dx::Mode::Tuned;
    else if (mode == "sig5") o.mode = lex1dx::Mode::Sig5;
    else throw std::invalid_argument("mode must be 'tuned' or 'sig5'");
    o.coarse_groups = M;
    lex1dx::Result r;
    {
        nb::gil_scoped_release release;
        r = lex1dx::cluster(x.data(), x.shape(0), k, o);
    }
    return to_dict(r);
}

}  // namespace

NB_MODULE(_core, m) {
    m.doc() = "lex1dx C++ core (internal; use the lex1dx package)";
    m.def("cluster", &cluster_fit, nb::arg("x"), nb::arg("k"), nb::arg("mode") = "tuned",
          nb::arg("M") = 0, nb::arg("tolerance") = 1e-11);
    m.def("default_M", &lex1dx::default_M);
}
