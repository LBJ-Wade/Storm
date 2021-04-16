#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <storm/simple.hpp>
#include <string>

namespace py = pybind11;

PYBIND11_MODULE(_simple_spectra, m) {
  m.doc() = "Module containing the spectra for the `SimpleRhNeutrino` model.";

  //================================
  //---- Two-Body Decay Spectra ----
  //================================

  m.def("spectrum_vr_to_vl_h", &storm::spectrum_vr_to_vl_h,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into an active neutrino and Higgs",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("product"),
        py::arg("xbounds"), py::arg("nbins"), py::arg("nevents"));

  m.def("spectrum_vr_to_vl_z", &storm::spectrum_vr_to_vl_z,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into an active neutrino and Z",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("product"),
        py::arg("xbounds"), py::arg("nbins"), py::arg("nevents"));

  m.def("spectrum_vr_to_l_w", &storm::spectrum_vr_to_l_w,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into a charged lepton and W^+",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("product"),
        py::arg("xbounds"), py::arg("nbins"), py::arg("nevents"),
        py::arg("anti"));

  //==================================
  //---- Three-Body Decay Spectra ----
  //==================================

  m.def("spectrum_vr_to_vl_u_u", &storm::spectrum_vr_to_vl_u_u,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into an active neutrino and two up-type quarks.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genq"),
        py::arg("product"), py::arg("xbounds"), py::arg("nbins"),
        py::arg("nevents"));

  m.def("spectrum_vr_to_vl_d_d", &storm::spectrum_vr_to_vl_d_d,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into an active neutrino and two down-type quarks.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genq"),
        py::arg("product"), py::arg("xbounds"), py::arg("nbins"),
        py::arg("nevents"));

  m.def("spectrum_vr_to_l_u_d", &storm::spectrum_vr_to_l_u_d,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into a charged lepton, an up-type quark and an anti-down-type quark.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genq"),
        py::arg("product"), py::arg("xbounds"), py::arg("nbins"),
        py::arg("nevents"), py::arg("anti"));

  m.def("spectrum_vr_to_vl_l_l", &storm::spectrum_vr_to_vl_l_l,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into an active neutrino and two charged leptons, all of the same "
        "generation",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("product"),
        py::arg("xbounds"), py::arg("nbins"), py::arg("nevents"));

  m.def("spectrum_vr_to_vl_lp_lp", &storm::spectrum_vr_to_vl_lp_lp,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into an active neutrino and two charged leptons of a different "
        "generation than neutrinos.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genlp"),
        py::arg("product"), py::arg("xbounds"), py::arg("nbins"),
        py::arg("nevents"));

  m.def("spectrum_vr_to_vlp_lp_l", &storm::spectrum_vr_to_vlp_lp_l,
        "Compute the spectrum of a given product from a RH neutrino decaying "
        "into an active neutrino and two charged lepton. Here the active "
        "neutrino and one of the leptons are of a different genertion than RH "
        "neutrino and one of the leptons is of the same generation as RH "
        "neutrino.",
        py::arg("mvr"), py::arg("theta"), py::arg("genl"), py::arg("genlp"),
        py::arg("product"), py::arg("xbounds"), py::arg("nbins"),
        py::arg("nevents"), py::arg("anti"));
}