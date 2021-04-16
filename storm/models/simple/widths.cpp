#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <storm/simple.hpp>

namespace py = pybind11;

PYBIND11_MODULE(_simple_widths, m) {
  m.doc() =
      "Module containing the partial widths for the `SimpleRhNeutrino` model.";

  m.def("width_vr_to_vl_h", &storm::width_vr_to_vl_h,
        "Compute the partial width for a RH neutrino decaying into an active "
        "neutrino and Higgs.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"));

  m.def("width_vr_to_vl_z", &storm::width_vr_to_vl_z,
        "Compute the partial width for a RH neutrino decaying into an active "
        "neutrino and Z.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"));

  m.def("width_vr_to_l_w", &storm::width_vr_to_l_w,
        "Compute the partial width for a RH neutrino decaying into a charged "
        "lepton and W.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"));

  m.def("width_vr_to_vl_u_u", &storm::width_vr_to_vl_u_u,
        "Compute the partial width for a RH neutrino decaying into an active "
        "neutrino and two up-type quarks.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genq"));

  m.def("width_vr_to_vl_d_d", &storm::width_vr_to_vl_d_d,
        "Compute the partial width for a RH neutrino decaying into an active "
        "neutrino and two down-type quarks.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genq"));

  m.def("width_vr_to_l_u_d", &storm::width_vr_to_l_u_d,
        "Compute the partial width for a RH neutrino decaying into an charged "
        "lepton, an up-type quark and a down-type quark.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genq"));

  m.def("width_vr_to_vl_lp_lp", &storm::width_vr_to_vl_lp_lp,
        "Compute the partial width for a RH neutrino decaying into an active "
        "neutrino and two charged leptons of a different generation that "
        "neutrinos.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genlp"));

  m.def("width_vr_to_vlp_lp_l", &storm::width_vr_to_vlp_lp_l,
        "Compute the partial decay width for a right-handed neutrino to "
        "decay into an active neutrino and two charged leptons, one of the "
        "same generation as the neutrino and the other a different generation.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genlp"));

  m.def("width_vr_to_vl_l_l", &storm::width_vr_to_vl_l_l,
        "Compute the partial decay width for a right-handed neutrino to "
        "decay into an active neutrino and two charged leptons, all of the "
        "same generation.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"));

  m.def("width_vr_to_vl_vl_vl", &storm::width_vr_to_vl_vl_vl,
        "Compute the partial decay width for a right-handed neutrino to "
        "decay into a three active neutrinos, all of the same generation.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"));
}