#include "storm/constants.hpp"
#include "storm/pythia_spectra.hpp"
#include "storm/simple.hpp"
#include <Pythia8/Basics.h>
#include <Pythia8/Pythia.h>
#include <functional>
#include <utility>
#include <vector>

namespace storm {

/**
 * Returns the PDG code of the neutrino given its generation.
 * @param gen Generation of the neutrino. Should be 0, 1 or 2.
 */
inline static auto pdg_nu(int gen) -> int { return 12 + 2 * gen; }

/**
 * Returns the PDG code of the charged lepton given its generation.
 * @param gen Generation of the lepton. Should be 0, 1 or 2.
 */
inline static auto pdg_lep(int gen) -> int { return 11 + 2 * gen; }

/**
 * Returns the PDG code of the up-type quark given its generation.
 * @param gen Generation of the up-type quark. Should be 0, 1 or 2.
 */
inline static auto pdg_qu(int gen) -> int { return 2 + 2 * gen; }

/**
 * Returns the PDG code of the up-type quark given its generation.
 * @param gen Generation of the up-type quark. Should be 0, 1 or 2.
 */
inline static auto pdg_qd(int gen) -> int { return 1 + 2 * gen; }

//=================================
//----- Two-Body Decay Spectra ----
//=================================

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and Higgs then to a specified final state.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 */
auto spectrum_vr_to_vl_h(double mvr, double theta, int gen, int product,
                         const std::pair<double, double> &xbounds,
                         unsigned int nbins, size_t nevents)
    -> std::pair<std::vector<double>, std::vector<double>> {

  if (mvr < kHIGGS_MASS) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }

  auto msqrd = [mvr, theta](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_vl_h(momenta, mvr, theta);
  };
  std::vector<int> ids{pdg_nu(gen), kHIGGS_ID};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and Z then to a specified final state.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 */
auto spectrum_vr_to_vl_z(double mvr, double theta, int gen, int product,
                         const std::pair<double, double> &xbounds,
                         unsigned int nbins, size_t nevents)
    -> std::pair<std::vector<double>, std::vector<double>> {

  if (mvr < kZ_BOSON_MASS) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }

  auto msqrd = [mvr, theta](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_vl_z(momenta, mvr, theta);
  };
  std::vector<int> ids{pdg_nu(gen), kZ_BOSON_ID};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

/**
 * Compute the spectrum from the decay of a right-handed neutrino into a charged
 * lepton and W^+ then to a specified final state.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param gen Generation of the final state lepton. Should be 0, 1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 * @param anti If true, use lbar + wbar.
 */
auto spectrum_vr_to_l_w(double mvr, double theta, int gen, int product,
                        const std::pair<double, double> &xbounds,
                        unsigned int nbins, size_t nevents, bool anti)
    -> std::pair<std::vector<double>, std::vector<double>> {

  const double ml = kLEPTON_MASSES.at(gen);
  if (mvr < ml + kW_BOSON_MASS) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }

  auto msqrd = [mvr, theta, gen](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_l_w(momenta, mvr, theta, gen);
  };
  const int eta = anti ? -1. : 1.0;
  std::vector<int> ids{eta * pdg_lep(gen), eta * kW_BOSON_ID};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

//===================================
//----- Three-Body Decay Spectra ----
//===================================

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and two up-type quarks.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 * @param genq Generation of the final state quarks. Should be 0, 1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 */
auto spectrum_vr_to_vl_u_u(double mvr, double theta, int genl, int genq,
                           int product,
                           const std::pair<double, double> &xbounds,
                           unsigned int nbins, size_t nevents)
    -> std::pair<std::vector<double>, std::vector<double>> {

  const double mu = kUP_QUARK_MASSES.at(genq);
  if (mvr < 2.0 * mu) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }

  auto msqrd = [mvr, theta, genq](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_vl_u_u(momenta, mvr, theta, genq);
  };
  std::vector<int> ids{pdg_nu(genl), pdg_qu(genq), -pdg_qu(genq)};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and two down-type quarks.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 * @param genq Generation of the final state quarks. Should be 0, 1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 */
auto spectrum_vr_to_vl_d_d(double mvr, double theta, int genl, int genq,
                           int product,
                           const std::pair<double, double> &xbounds,
                           unsigned int nbins, size_t nevents)
    -> std::pair<std::vector<double>, std::vector<double>> {

  const double md = kDOWN_QUARK_MASSES.at(genq);
  if (mvr < 2.0 * md) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }

  auto msqrd = [mvr, theta, genq](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_vl_d_d(momenta, mvr, theta, genq);
  };
  std::vector<int> ids{pdg_nu(genl), pdg_qd(genq), -pdg_qd(genq)};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

/**
 * Compute the spectrum from the decay of a right-handed neutrino into a charged
 * lepton, an up-type quark and a anti-down-type quark.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 * @param genq Generation of the final state quarks. Should be 0, 1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 * @param anti If true, use lbar + ubar + d.
 */
auto spectrum_vr_to_l_u_d(double mvr, double theta, int genl, int genq,
                          int product, const std::pair<double, double> &xbounds,
                          unsigned int nbins, size_t nevents, bool anti)
    -> std::pair<std::vector<double>, std::vector<double>> {

  const double ml = kLEPTON_MASSES.at(genl);
  const double mu = kUP_QUARK_MASSES.at(genq);
  const double md = kDOWN_QUARK_MASSES.at(genq);
  if (mvr < ml + mu + md) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }

  auto msqrd = [mvr, theta, genl,
                genq](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_l_u_d(momenta, mvr, theta, genl, genq);
  };
  const int eta = anti ? -1.0 : 1.0;
  std::vector<int> ids{eta * pdg_lep(genl), eta * pdg_qu(genq),
                       -eta * pdg_qd(genq)};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and two charged leptons, all of the same generation.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 */
auto spectrum_vr_to_vl_l_l(double mvr, double theta, int genl, int product,
                           const std::pair<double, double> &xbounds,
                           unsigned int nbins, size_t nevents)
    -> std::pair<std::vector<double>, std::vector<double>> {

  const double ml = kLEPTON_MASSES.at(genl);
  if (mvr < 2 * ml) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }

  auto msqrd = [mvr, theta, genl](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_vl_l_l(momenta, mvr, theta, genl);
  };
  std::vector<int> ids{pdg_nu(genl), pdg_lep(genl), -pdg_lep(genl)};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and two charged leptons, where the charge leptons are of a different
 * generation than active neutrino.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 * @param genlp Generation of the final state chaged lepton. Should be 0,1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 */
auto spectrum_vr_to_vl_lp_lp(double mvr, double theta, int genl, int genlp,
                             int product,
                             const std::pair<double, double> &xbounds,
                             unsigned int nbins, size_t nevents)
    -> std::pair<std::vector<double>, std::vector<double>> {

  const double mlp = kLEPTON_MASSES.at(genlp);
  if (mvr < 2 * mlp) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }

  auto msqrd = [mvr, theta, genl,
                genlp](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_vl_lp_lp(momenta, mvr, theta, genl, genlp);
  };
  std::vector<int> ids{pdg_nu(genl), pdg_lep(genlp), -pdg_lep(genlp)};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and two charged leptons, where the charge leptons are of a different
 * generation than active neutrino.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 * @param genlp Generation of the final state chaged lepton. Should be 0,1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 * @param anti If true, use vlp + lpbar + lp.
 */
auto spectrum_vr_to_vlp_lp_l(double mvr, double theta, int genl, int genlp,
                             int product,
                             const std::pair<double, double> &xbounds,
                             unsigned int nbins, size_t nevents, bool anti)
    -> std::pair<std::vector<double>, std::vector<double>> {

  const double ml = kLEPTON_MASSES.at(genl);
  const double mlp = kLEPTON_MASSES.at(genlp);
  if (mvr < ml + mlp) {
    return std::make_pair(std::vector<double>{}, std::vector<double>{});
  }
  auto msqrd = [mvr, theta, genl,
                genlp](const std::vector<Pythia8::Vec4> &momenta) {
    return msqrd_vr_to_vlp_lp_l(momenta, mvr, theta, genl, genlp);
  };
  const int eta = anti ? -1.0 : 1.0;
  std::vector<int> ids{pdg_nu(genl), eta * pdg_lep(genlp),
                       -eta * pdg_lep(genl)};
  return decay_spectrum(xbounds, nbins, mvr, product, ids, msqrd, nevents);
}

} // namespace storm