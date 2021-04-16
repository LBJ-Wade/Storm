/**
 * Header file containing all the definitions for the simple right-handed (RH)
 * neutrino model where the RH neutrino couples to a single active neutrino
 * through a Yukawa interaction with Higgs.
 */

#ifndef STORM_SIMPLE_HPP
#define STORM_SIMPLE_HPP

#include "storm/types.hpp"
#include <Pythia8/Basics.h>
#include <Pythia8/Pythia.h>
#include <cstddef>
#include <functional>
#include <locale>
#include <map>
#include <vector>

namespace storm {

// ========================
// ---- Partial Widths ----
// ========================

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino and Higgs.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 */
auto width_vr_to_vl_h(double, double, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino and Z.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 */
auto width_vr_to_vl_z(double, double, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into a charged
 * lepton and W.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 */
auto width_vr_to_l_w(double, double, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino and two up-type quarks.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the lepton.
 * @param genq Generation of the up-type quarks.
 * @returns The width.
 */
auto width_vr_to_vl_u_u(double, double, int, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino and two down-type quarks.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the lepton.
 * @param genq Generation of the down-type quarks.
 * @returns The width.
 */
auto width_vr_to_vl_d_d(double, double, int, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into a charged
 * lepton, an up-type quark and a down-type quark.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the lepton.
 * @param genq Generation of the quarks.
 * @returns The width.
 */
auto width_vr_to_l_u_d(double, double, int, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino, and two charged leptons with different generations than neutrino.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the neutrino.
 * @param genlp Generation of charged leptons.
 * @returns The width.
 */
auto width_vr_to_vl_lp_lp(double, double, int, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino, and two charged leptons with different generations than neutrino.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the neutrino.
 * @param genlp Generation of charged leptons.
 * @returns The width.
 */
auto width_vr_to_vlp_lp_l(double, double, int, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino, and two charged leptons with of the same generation as neutrino.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the final-state leptons.
 * @returns The width.
 */
auto width_vr_to_vl_l_l(double, double, int) -> double;

/**
 * Compute the partial width for a right-handed neutrino to decay into a three
 * active neutrinos, all of the same generation.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the final-state leptons.
 * @returns The width.
 */
auto width_vr_to_vl_vl_vl(double, double, int) -> double;

// =================================
// ---- Squared Matrix Elements ----
// =================================

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and Higgs.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 */
auto msqrd_vr_to_vl_h(const std::vector<Pythia8::Vec4> &, double, double)
    -> double;

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and Z.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 */
auto msqrd_vr_to_vl_z(const std::vector<Pythia8::Vec4> &, double, double)
    -> double;

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * a charged lepton and W.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the final state lepton. Should be 0, 1 or 2.
 */
auto msqrd_vr_to_l_w(const std::vector<Pythia8::Vec4> &, double, double, int)
    -> double;

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and two up-type quarks.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genq Generation of the up-type quarks. Should be 0, 1 or 2.
 */
auto msqrd_vr_to_vl_u_u(const std::vector<Pythia8::Vec4> &, double, double, int)
    -> double;

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and two down-type quarks.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genq Generation of the up-type quarks. Should be 0, 1 or 2.
 */
auto msqrd_vr_to_vl_d_d(const std::vector<Pythia8::Vec4> &, double, double, int)
    -> double;

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * a charged lepton, an up-type quark and a down-type quark.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the charged lepton.
 * @param genq Generation of the quarks.
 */
auto msqrd_vr_to_l_u_d(const std::vector<Pythia8::Vec4> &, double, double, int,
                       int) -> double;

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and two charged leptons of a different generation than
 * neutrino.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the neutrinos.
 * @param genlp Generation of the charged leptons.
 */
auto msqrd_vr_to_vl_lp_lp(const std::vector<Pythia8::Vec4> &, double, double,
                          int, int) -> double;

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and two charged leptons, one of the same generation as the
 * neutrino and the other a different generation.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the right-handed neutrino.
 * @param genlp Generation of the other leptons.
 */
auto msqrd_vr_to_vlp_lp_l(const std::vector<Pythia8::Vec4> &, double, double,
                          int, int) -> double;

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and two charged leptons, all of the same generation.
 * neutrino.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the leptons.
 */
auto msqrd_vr_to_vl_l_l(const std::vector<Pythia8::Vec4> &, double, double, int)
    -> double;

// ========================================
// ---- Partially Integraded Functions ----
// ========================================

struct PartiallyIntegratedMsqrdParams {
  const double mvr;
  const double theta;
  const int genl;
  const int genq;
};

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and up-type quarks.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vl_u_u(double, void *) -> double;

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and down-type quarks.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vl_d_d(double, void *) -> double;

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into a charged lepton, an up-type quark and a down-type
 * quark.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_l_u_d(double, void *) -> double;

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and two charged leptons of
 * different generation than neutrinos.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vl_lp_lp(double, void *) -> double;

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and two charged leptons. Here
 * the neutrino shares generation with the charged lepton and the right-haned
 * neutrino is of the same generation as anti-charged lepton.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vlp_lp_l(double, void *) -> double;

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and two charged lepton, all
 * with the same generation. quark.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vl_l_l(double, void *) -> double;

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an three active neutrino, all of the same
 * generation.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta.
 */
auto partialy_integrated_msqrd_vr_to_vl_vl_vl(double, void *) -> double;

// =======================
// ---- Decay Spectra ----
// =======================

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and Higgs then to a specified final state.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param gen Generation of the final state lepton. Should be 0, 1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 */
auto spectrum_vr_to_vl_h(double, double, int, int,
                         const std::pair<double, double> &, unsigned int,
                         size_t)
    -> std::pair<std::vector<double>, std::vector<double>>;

/**
 * Compute the spectrum from the decay of a right-handed neutrino into an active
 * neutrino and Z then to a specified final state.
 * @param mvr Mass of the right-handed neutrino.
 * @param theta Mixing angle between right-handed neutrino and active neutrino.
 * @param gen Generation of the final state lepton. Should be 0, 1 or 2.
 * @param product PDG code of the product to compute spectrum for.
 * @param xbounds Bounds on x = 2*E/mvr.
 * @param nbins Number of bins/energies to compute spectrum for.
 * @param nevents Number of pythia events to run.
 */
auto spectrum_vr_to_vl_z(double, double, int, int,
                         const std::pair<double, double> &, unsigned int,
                         size_t)
    -> std::pair<std::vector<double>, std::vector<double>>;

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
auto spectrum_vr_to_l_w(double, double, int, int,
                        const std::pair<double, double> &, unsigned int, size_t,
                        bool)
    -> std::pair<std::vector<double>, std::vector<double>>;

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
auto spectrum_vr_to_vl_u_u(double, double, int, int, int,
                           const std::pair<double, double> &, unsigned int,
                           size_t)
    -> std::pair<std::vector<double>, std::vector<double>>;

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
auto spectrum_vr_to_vl_d_d(double, double, int, int, int,
                           const std::pair<double, double> &, unsigned int,
                           size_t)
    -> std::pair<std::vector<double>, std::vector<double>>;

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
auto spectrum_vr_to_l_u_d(double, double, int, int, int,
                          const std::pair<double, double> &, unsigned int,
                          size_t, bool)
    -> std::pair<std::vector<double>, std::vector<double>>;

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
auto spectrum_vr_to_vl_l_l(double, double, int, int,
                           const std::pair<double, double> &, unsigned int,
                           size_t)
    -> std::pair<std::vector<double>, std::vector<double>>;

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
auto spectrum_vr_to_vl_lp_lp(double, double, int, int, int,
                             const std::pair<double, double> &, unsigned int,
                             size_t)
    -> std::pair<std::vector<double>, std::vector<double>>;

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
auto spectrum_vr_to_vlp_lp_l(double, double, int, int, int,
                             const std::pair<double, double> &, unsigned int,
                             size_t, bool)
    -> std::pair<std::vector<double>, std::vector<double>>;

} // namespace storm

#endif // STORM_SIMPLE_HPP
