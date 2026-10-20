/**
 * @file bindings.cpp
 * @brief pybind11 module `lsquant._core` exposing the C++ Chebyshev engine to Python.
 *
 * Binds `SparseMatrix`, the `LinearOperator` hierarchy, `StateGenerator` implementations,
 * `SpectralWindow`, `RescaledOperator`, the five bytecode instructions and `Plan`,
 * `MomentEngine`, `ChebyshevMoments` (with a zero-copy NumPy `complex128` view), file I/O
 * routines, and rank-1 series contraction (`ContractRank1`, `WriteCurveText`).
 */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "core/algebra/composite_operators.hpp"
#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/state_generator.hpp"
#include "core/algebra/types.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"
#include "core/chebyshev/contraction.hpp"
#include "core/chebyshev/moment_engine.hpp"
#include "core/chebyshev/plan.hpp"
#include "core/chebyshev/rescaled_operator.hpp"
#include "core/chebyshev/spectral_window.hpp"
#include "io/chebmom_file.hpp"
#include "io/csr_text.hpp"
#include "io/curve_text.hpp"

namespace py = pybind11;

namespace lsquant {
namespace {

/**
 * @brief Validate a 2D `float64` NumPy array and contract it with a strided real moment vector.
 *
 * @param momentsPtr   Pointer to the first real moment element.
 * @param momentStride Element stride between consecutive moments.
 * @param numMoments   Number of moments `M`.
 * @param table        2D `float64` NumPy array of shape `(M, N)`.
 * @return             1D `float64` NumPy array of length `N`.
 */
py::array_t<double> ContractRank1Buffer(const double* momentsPtr, Index momentStride,
                                        Index numMoments, const py::array_t<double>& table) {
  if (table.ndim() != 2) {
    throw std::invalid_argument("ContractRank1: table must be a 2D float64 array of shape (M, N)");
  }
  if (static_cast<Index>(table.shape(0)) != numMoments) {
    throw std::invalid_argument(
        "ContractRank1: table.shape[0] must equal the number of moments M");
  }
  const Index numPoints = static_cast<Index>(table.shape(1));
  const auto elemSize = static_cast<py::ssize_t>(sizeof(double));
  if (table.strides(0) % elemSize != 0 || table.strides(1) % elemSize != 0) {
    throw std::invalid_argument("ContractRank1: unaligned table strides");
  }
  const Index tableStrideM = static_cast<Index>(table.strides(0) / elemSize);
  const Index tableStrideN = static_cast<Index>(table.strides(1) / elemSize);

  py::array_t<double> result(static_cast<py::ssize_t>(numPoints));
  ContractRank1(momentsPtr, momentStride, numMoments, table.data(), tableStrideM, tableStrideN,
                numPoints, result.mutable_data());
  return result;
}

}  // namespace
}  // namespace lsquant

PYBIND11_MODULE(_core, m) {
  using namespace lsquant;

  m.doc() = "C++ execution engine and zero-copy tensor bindings for LSQUANT.";

  py::class_<SparseMatrix>(m, "SparseMatrix")
      .def(py::init<>())
      .def_static(
          "from_csr",
          [](Index dim, const py::array_t<Index, py::array::c_style>& rowPointers,
             const py::array_t<Index, py::array::c_style>& columns,
             const py::array_t<Scalar, py::array::c_style>& values) {
            if (rowPointers.ndim() != 1 || columns.ndim() != 1 || values.ndim() != 1) {
              throw std::invalid_argument(
                  "SparseMatrix.from_csr: row_pointers, columns, and values must be 1D arrays");
            }
            std::vector<Index> rp(rowPointers.data(), rowPointers.data() + rowPointers.size());
            std::vector<Index> cols(columns.data(), columns.data() + columns.size());
            std::vector<Scalar> vals(values.data(), values.data() + values.size());
            return SparseMatrix::FromCsr(dim, std::move(rp), std::move(cols), std::move(vals));
          },
          py::arg("dim"), py::arg("row_pointers"), py::arg("columns"), py::arg("values"))
      .def_property_readonly("dimension", &SparseMatrix::Dimension)
      .def_property_readonly("non_zeros", &SparseMatrix::NonZeros)
      .def("hermitian_fold", &SparseMatrix::HermitianFold)
      .def("is_hermitian", &SparseMatrix::IsHermitian, py::arg("tol") = 1e-12);

  py::class_<LinearOperator>(m, "LinearOperator")
      .def_property_readonly("dimension", &LinearOperator::Dimension)
      .def_property_readonly("is_identity", &LinearOperator::IsIdentity)
      .def_property_readonly("is_hermitian", &LinearOperator::IsHermitian);

  py::class_<SparseMatrixOperator, LinearOperator>(m, "SparseMatrixOperator")
      .def(py::init<const SparseMatrix&>(), py::arg("matrix"), py::keep_alive<1, 2>());

  py::class_<IdentityOperator, LinearOperator>(m, "IdentityOperator")
      .def(py::init<Index>(), py::arg("dim"));

  py::class_<StateGenerator>(m, "StateGenerator")
      .def_property_readonly("dimension", &StateGenerator::Dimension)
      .def_property_readonly("number_of_states", &StateGenerator::NumberOfStates)
      .def_property_readonly("label", &StateGenerator::Label)
      .def("reset", &StateGenerator::Reset);

  py::class_<RandomPhaseStates, StateGenerator>(m, "RandomPhaseStates")
      .def(py::init<Index, Index, unsigned>(), py::arg("dim"), py::arg("num_states"),
           py::arg("seed") = 12345u);

  py::class_<LocalStates, StateGenerator>(m, "LocalStates")
      .def(py::init([](Index dim, const py::array_t<Index, py::array::c_style>& sites) {
             if (sites.ndim() != 1) {
               throw std::invalid_argument("LocalStates: sites must be a 1D int64 array");
             }
             std::vector<Index> siteVec(sites.data(), sites.data() + sites.size());
             return LocalStates(dim, std::move(siteVec));
           }),
           py::arg("dim"), py::arg("sites"))
      .def(py::init<Index, std::vector<Index>>(), py::arg("dim"), py::arg("sites"));

  py::class_<UserStates, StateGenerator>(m, "UserStates")
      .def(py::init([](const std::vector<py::array_t<Scalar, py::array::c_style>>& arrays) {
             std::vector<Vector> states;
             states.reserve(arrays.size());
             for (const auto& arr : arrays) {
               if (arr.ndim() != 1) {
                 throw std::invalid_argument("UserStates: each state must be a 1D complex128 array");
               }
               states.emplace_back(arr.data(), arr.data() + arr.size());
             }
             return UserStates(std::move(states));
           }),
           py::arg("states"))
      .def(py::init<std::vector<Vector>>(), py::arg("states"));

  py::class_<SpectralWindow>(m, "SpectralWindow")
      .def(py::init<>())
      .def(py::init<double, double>(), py::arg("lower"), py::arg("upper"))
      .def_readwrite("lower", &SpectralWindow::lower)
      .def_readwrite("upper", &SpectralWindow::upper)
      .def_property_readonly("scale", &SpectralWindow::Scale)
      .def_property_readonly("shift", &SpectralWindow::Shift)
      .def_property_readonly("half_width", &SpectralWindow::HalfWidth)
      .def_property_readonly("center", &SpectralWindow::Center);

  py::class_<RescaledOperator, LinearOperator>(m, "RescaledOperator")
      .def(py::init<const SparseMatrix&, const SpectralWindow&>(), py::arg("matrix"),
           py::arg("window"), py::keep_alive<1, 2>())
      .def_property_readonly("window", &RescaledOperator::Window);

  py::class_<ApplyOperator>(m, "ApplyOperator")
      .def(py::init<int>(), py::arg("operator_index"))
      .def_readwrite("operator_index", &ApplyOperator::operatorIndex);

  py::class_<Expand>(m, "Expand")
      .def(py::init<int>(), py::arg("num_moments"))
      .def_readwrite("num_moments", &Expand::numMoments);

  py::class_<ApplySeries>(m, "ApplySeries")
      .def(py::init([](const py::array_t<Scalar, py::array::c_style>& coeffs) {
             if (coeffs.ndim() != 1) {
               throw std::invalid_argument("ApplySeries: coeffs must be a 1D complex128 array");
             }
             return ApplySeries{std::vector<Scalar>(coeffs.data(), coeffs.data() + coeffs.size())};
           }),
           py::arg("coeffs"))
      .def(py::init<std::vector<Scalar>>(), py::arg("coeffs"))
      .def_readwrite("coeffs", &ApplySeries::coeffs);

  py::class_<MirrorToBra>(m, "MirrorToBra").def(py::init<>());

  py::class_<Trace>(m, "Trace")
      .def(py::init([](const py::array_t<Scalar, py::array::c_style>& weights) {
             if (weights.ndim() != 1) {
               throw std::invalid_argument("Trace: weights must be a 1D complex128 array");
             }
             return Trace{std::vector<Scalar>(weights.data(), weights.data() + weights.size())};
           }),
           py::arg("weights"))
      .def(py::init<std::vector<Scalar>>(), py::arg("weights"))
      .def_readwrite("weights", &Trace::weights);

  py::class_<Plan>(m, "Plan")
      .def(py::init<>())
      .def(py::init<std::vector<Instruction>, std::string>(),
           py::arg("code") = std::vector<Instruction>{}, py::arg("word") = "")
      .def_readwrite("code", &Plan::code)
      .def_readwrite("word", &Plan::word);

  py::class_<MomentEngine>(m, "MomentEngine")
      .def(py::init<const LinearOperator&, std::vector<const LinearOperator*>, StateGenerator&,
                    int>(),
           py::arg("rescaled_h"), py::arg("operators"), py::arg("states"), py::arg("batch_size"),
           py::keep_alive<1, 2>(), py::keep_alive<1, 3>(), py::keep_alive<1, 4>())
      .def("run", &MomentEngine::Run, py::arg("plan"));

  py::class_<ChebyshevMoments>(m, "ChebyshevMoments")
      .def(py::init<std::vector<int>>(), py::arg("dims") = std::vector<int>{})
      .def_property_readonly("rank", &ChebyshevMoments::Rank)
      .def_property_readonly("dims", &ChebyshevMoments::Dims)
      .def_property_readonly(
          "data",
          [](py::object selfObj) {
            auto& self = selfObj.cast<ChebyshevMoments&>();
            std::vector<py::ssize_t> shape;
            shape.reserve(self.Dims().size());
            for (int d : self.Dims()) {
              shape.push_back(static_cast<py::ssize_t>(d));
            }
            std::vector<py::ssize_t> strides(shape.size());
            py::ssize_t stride = static_cast<py::ssize_t>(sizeof(Scalar));
            for (py::ssize_t i = static_cast<py::ssize_t>(shape.size()) - 1; i >= 0; --i) {
              strides[static_cast<std::size_t>(i)] = stride;
              stride *= shape[static_cast<std::size_t>(i)];
            }
            return py::array_t<Scalar>(shape, strides, self.Data().data(), selfObj);
          })
      .def_readwrite("window", &ChebyshevMoments::window)
      .def_readwrite("system_size", &ChebyshevMoments::systemSize)
      .def_readwrite("number_of_states", &ChebyshevMoments::numberOfStates)
      .def_readwrite("state_label", &ChebyshevMoments::stateLabel)
      .def_readwrite("word", &ChebyshevMoments::word)
      .def_readwrite("kpm_weighted", &ChebyshevMoments::kpmWeighted)
      .def_readwrite("symmetrised", &ChebyshevMoments::symmetrised);

  m.def("ApplyKpmWeights", &ApplyKpmWeights, py::arg("moments"), py::arg("num_states"),
        py::arg("symmetrise"));

  m.def("ReadCsrText", &ReadCsrText, py::arg("path"));
  m.def("ReadChebmom", &ReadChebmom, py::arg("path"));
  m.def("WriteChebmom", &WriteChebmom, py::arg("path"), py::arg("moments"));
  m.def("WriteLegacyChebmom1D", &WriteLegacyChebmom1D, py::arg("path"), py::arg("moments"));
  m.def("WriteLegacyChebmom2D", &WriteLegacyChebmom2D, py::arg("path"), py::arg("moments"));

  m.def(
      "ContractRank1",
      [](const py::array_t<double>& moments, const py::array_t<double>& table) {
        if (moments.ndim() != 1) {
          throw std::invalid_argument("ContractRank1: moments must be a 1D float64 array");
        }
        const auto elemSize = static_cast<py::ssize_t>(sizeof(double));
        if (moments.strides(0) % elemSize != 0) {
          throw std::invalid_argument("ContractRank1: unaligned moments stride");
        }
        const Index momentStride = static_cast<Index>(moments.strides(0) / elemSize);
        const Index numMoments = static_cast<Index>(moments.shape(0));
        return ContractRank1Buffer(moments.data(), momentStride, numMoments, table);
      },
      py::arg("moments"), py::arg("table"));

  m.def(
      "ContractRank1",
      [](const ChebyshevMoments& moments, const py::array_t<double>& table) {
        if (moments.Rank() != 1) {
          throw std::invalid_argument("ContractRank1: ChebyshevMoments must have rank 1");
        }
        const Index numMoments = static_cast<Index>(moments.Dims()[0]);
        const auto* realPtr = reinterpret_cast<const double*>(moments.Data().data());
        return ContractRank1Buffer(realPtr, 2, numMoments, table);
      },
      py::arg("moments"), py::arg("table"));

  m.def(
      "WriteCurveText",
      [](const std::string& path, const py::array_t<double>& x, const py::array_t<double>& y) {
        if (x.ndim() != 1 || y.ndim() != 1) {
          throw std::invalid_argument("WriteCurveText: x and y must be 1D float64 arrays");
        }
        if (x.shape(0) != y.shape(0)) {
          throw std::invalid_argument("WriteCurveText: x and y must have the same length");
        }
        const auto elemSize = static_cast<py::ssize_t>(sizeof(double));
        if (x.strides(0) % elemSize != 0 || y.strides(0) % elemSize != 0) {
          throw std::invalid_argument("WriteCurveText: unaligned array strides");
        }
        const Index xStride = static_cast<Index>(x.strides(0) / elemSize);
        const Index yStride = static_cast<Index>(y.strides(0) / elemSize);
        WriteCurveText(path, x.data(), xStride, y.data(), yStride,
                       static_cast<Index>(x.shape(0)));
      },
      py::arg("path"), py::arg("x"), py::arg("y"));
}
