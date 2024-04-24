#include <pybind11/pybind11.h>
#include <matplot/matplot.h>

namespace py = pybind11;
namespace mp = matplot;

float square(float x) { return x * x; }

PYBIND11_MODULE(projekt, m) {
    m.def("square", &square);
}