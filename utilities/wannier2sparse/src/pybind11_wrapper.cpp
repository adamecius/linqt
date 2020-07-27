
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


#include<vector>
#include "wannier2sparse.hpp"
using namespace std;
int py_wannier2sparse(string label, array<int, 3> cellDim, deque< string > op_list)
{
	return wannier2sparse(label, cellDim, op_list);
}


namespace py = pybind11;

PYBIND11_MODULE(wannier2sparse, m) {
    m.doc() = R"pbdoc(
        Pybind11 wrapper of wannier2sparse
        -----------------------
        .. currentmodule:: wannier2sparse
        .. autosummary::
           :toctree: _generate
           wannier2sparse
    )pbdoc";

    m.def("wannier2sparse", &py_wannier2sparse, R"pbdoc(
		This function map a Wannier Hamiltonian in the Wannier90 format
		into a supercell in real, momentum or eigenvector spaces which is
		stored as an sparse matrix.
		)pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
