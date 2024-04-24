#include <pybind11/pybind11.h>
#include <matplot/matplot.h>

namespace py = pybind11;
namespace mp = matplot;

//testowa funkcja do usunięcia po prostu sprawdzałem czy się kompiluje jak należy
float square(float x) { return x * x; }
//////////////////////////////////////
PYBIND11_MODULE(projekt, m) {
    //również do usunięcia póki co zostawiam jako template
    m.def("square", &square);
    ///////////////////////////////////
}
