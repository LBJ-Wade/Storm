/**
 * This file implements the partial widths for the `SimpleRhNeutrino` model.
 */

#include "storm/constants.hpp"
#include "storm/rambo.hpp"
#include "storm/simple.hpp"
#include "storm/types.hpp"
#include <Pythia8/Basics.h>
#include <algorithm>
#include <array>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <cstddef>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <locale>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace storm {

//=========================
//---- Two-Body Widths ----
//=========================

auto width_vr_to_vl_h(double mvr, double theta, int /*genl*/) -> double {
  if (mvr < kHIGGS_MASS) {
    return 0.0;
  }
  return ((-pow(kHIGGS_MASS, 2) + pow(mvr, 2)) * pow(theta, 2) *
          std::abs(pow(kHIGGS_MASS, 2) - pow(mvr, 2))) /
         (16. * mvr * M_PI * pow(kHIGGS_VEV, 2));
}

auto width_vr_to_vl_z(double mvr, double theta, int /*genl*/) -> double {
  if (mvr < kZ_BOSON_MASS) {
    return 0.0;
  }
  return (((pow(mvr, 4) + pow(mvr, 2) * pow(kZ_BOSON_MASS, 2) -
            2 * pow(kZ_BOSON_MASS, 4)) *
           kALPHA_EM * pow(theta, 2) *
           std::abs(pow(mvr, 2) - pow(kZ_BOSON_MASS, 2))) /
          (16.0 * pow(kCOS_THETA_WEAK, 2) * pow(mvr, 3) *
           pow(kZ_BOSON_MASS, 2) * pow(kSIN_THETA_WEAK, 2)));
}

auto width_vr_to_l_w(double mvr, double theta, int genl) -> double {

  if (genl < 0 || genl > 2) {
    // throw std::invalid_argument(
    //    "Invalid generation passed to width_vr_l_w. Use 0, 1 or 2.");
    return 0.0;
  }

  const double ml = kLEPTON_MASSES.at(genl);

  if (mvr < ml + kW_BOSON_MASS) {
    return 0.0;
  }

  return (std::sqrt((ml - mvr - kW_BOSON_MASS) * (ml + mvr - kW_BOSON_MASS) *
                    (ml - mvr + kW_BOSON_MASS) * (ml + mvr + kW_BOSON_MASS)) *
          (pow(pow(ml, 2) - pow(mvr, 2), 2) +
           (pow(ml, 2) + pow(mvr, 2)) * pow(kW_BOSON_MASS, 2) -
           2 * pow(kW_BOSON_MASS, 4)) *
          kALPHA_EM * pow(theta, 2)) /
         (16. * pow(mvr, 3) * pow(kW_BOSON_MASS, 2) * pow(kSIN_THETA_WEAK, 2));
}

//===========================
//---- Three-Body Widths ----
//===========================

/**
 * Compute the integration bounds on the Mandelstam variable `s` for three-body
 * phase-space integration.
 * @param m Mass of the decaying particle.
 * @param m1 Mass of final state particle 1.
 * @param m2 Mass of final state particle 2.
 * @param m3 Mass of final state particle 3.
 */
inline static auto compute_s_bounds(double m, double m1, double m2, double m3)
    -> std::pair<double, double> {
  return std::make_pair(pow(m2 + m3, 2), pow(m - m1, 2));
}

/**
 * Compute the integration bounds on the mandelstam variable t.
 */
inline static auto compute_t_bounds(double s, double m, double m1, double m2,
                                    double m3) -> std::pair<double, double> {
  const double t1 = (pow(m1, 2) * pow(m2, 2) - pow(m1, 2) * pow(m3, 2) +
                     pow(m1, 2) * s + pow(m2, 2) * s + pow(m3, 2) * s -
                     pow(s, 2) + pow(m, 2) * (-pow(m2, 2) + pow(m3, 2) + s) -
                     sqrt(pow(m, 4) + pow(pow(m1, 2) - s, 2) -
                          2 * pow(m, 2) * (pow(m1, 2) + s)) *
                         sqrt(pow(m2, 4) + pow(pow(m3, 2) - s, 2) -
                              2 * pow(m2, 2) * (pow(m3, 2) + s))) /
                    (2. * s);
  const double t2 = (pow(m1, 2) * pow(m2, 2) - pow(m1, 2) * pow(m3, 2) +
                     pow(m1, 2) * s + pow(m2, 2) * s + pow(m3, 2) * s -
                     pow(s, 2) + pow(m, 2) * (-pow(m2, 2) + pow(m3, 2) + s) +
                     sqrt(pow(m, 4) + pow(pow(m1, 2) - s, 2) -
                          2 * pow(m, 2) * (pow(m1, 2) + s)) *
                         sqrt(pow(m2, 4) + pow(pow(m3, 2) - s, 2) -
                              2 * pow(m2, 2) * (pow(m3, 2) + s))) /
                    (2. * s);
  return (t1 < t2) ? std::make_pair(t1, t2) : std::make_pair(t2, t1);
}

/**
 * Compute a three-body decay width given the squared matrix element, mass of
 * the decaying particle and final-state particle masses using RAMBO.
 */
static auto compute_width_3body(const SquaredMatrixElement &msqrd,
                                const double m, std::vector<double> fsp_masses,
                                size_t num_events)
    -> std::pair<double, double> {

  if (m < std::accumulate(fsp_masses.begin(), fsp_masses.end(), 0.0)) {
    return std::make_pair(0.0, 0.0);
  }

  Rambo rambo{std::move(fsp_masses), m, msqrd};
  const auto width = rambo.compute_width(num_events);
  return width;
}

/**
 * Compute a three-body decay width given a function for the
 * partially-integrated (over Mandelstam variable `t`).
 */
static auto compute_width_3body_quad(gsl_function *integrand, const double m,
                                     std::array<double, 3> fsp_masses)
    -> std::pair<double, double> {

  if (m < std::accumulate(fsp_masses.begin(), fsp_masses.end(), 0.0)) {
    return std::make_pair(0.0, 0.0);
  }

  gsl_set_error_handler_off();

  auto [s_low, s_high] =
      compute_s_bounds(m, fsp_masses[0], fsp_masses[1], fsp_masses[2]);

  // Break-points for integration
  std::vector<double> bpts;
  bpts.reserve(5);
  bpts.push_back(s_low);
  // If MW^2, MZ^2 or MH^2 is in the integration interval, add them to
  // break-points
  const double mw2 = kW_BOSON_MASS * kW_BOSON_MASS;
  const double mz2 = kZ_BOSON_MASS * kZ_BOSON_MASS;
  const double mh2 = kHIGGS_MASS * kHIGGS_MASS;
  if (s_low < mw2 && mw2 < s_high) {
    bpts.push_back(mw2);
  }
  if (s_low < mz2 && mz2 < s_high) {
    bpts.push_back(mz2);
  }
  if (s_low < mh2 && mh2 < s_high) {
    bpts.push_back(mh2);
  }
  bpts.push_back(s_high);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  double integral = 0.0;
  double error = 0.0;
  gsl_integration_qagp(integrand, bpts.data(), bpts.size(), 0.0, 1e-7, 1000, w,
                       &integral, &error);

  const double pf = 1.0 / (32.0 * pow(2.0 * M_PI * m, 3));
  return std::make_pair(std::abs(integral * pf), error * pf);
}

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino and two up-type quarks.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the lepton.
 * @param genq Generation of the down-type quarks.
 * @returns The width and an error estimate.
 */
auto width_vr_to_vl_u_u(double mvr, double theta, int genl, int genq)
    -> double {
  const double mu = kUP_QUARK_MASSES.at(genq);
  // auto msqrd = [mvr, theta, genq](const std::vector<Pythia8::Vec4> &momenta)
  // {
  //  return msqrd_vr_to_vl_u_u(momenta, mvr, theta, genq);
  //};
  // auto msqrd = [mvr, theta, genl, genq](double s) {
  //  return partialy_integrated_msqrd_vr_to_vl_u_u(s, mvr, theta, genl, genq);
  //};

  PartiallyIntegratedMsqrdParams params{mvr, theta, genl, genq};

  gsl_function integrand;
  integrand.function = &partialy_integrated_msqrd_vr_to_vl_u_u;
  integrand.params = &params;

  return compute_width_3body_quad(&integrand, mvr, {0.0, mu, mu}).first;
}

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino and two down-type quarks.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the lepton.
 * @param genq Generation of the up-type quarks.
 * @returns The width and an error estimate.
 */
auto width_vr_to_vl_d_d(double mvr, double theta, int genl, int genq)
    -> double {
  const double md = kDOWN_QUARK_MASSES.at(genq);
  // auto msqrd = [mvr, theta, genq](const std::vector<Pythia8::Vec4> &momenta)
  // {
  //  return msqrd_vr_to_vl_d_d(momenta, mvr, theta, genq);
  //};
  // return compute_width_3body(msqrd, mvr, {0.0, md, md}, nevents);

  // auto msqrd = [mvr, theta, genl, genq](double s) {
  //  return partialy_integrated_msqrd_vr_to_vl_u_u(s, mvr, theta, genl, genq);
  //};

  PartiallyIntegratedMsqrdParams params{mvr, theta, genl, genq};

  gsl_function integrand;
  integrand.function = &partialy_integrated_msqrd_vr_to_vl_d_d;
  integrand.params = &params;

  return compute_width_3body_quad(&integrand, mvr, {0.0, md, md}).first;
}

/**
 * Compute the partial width for a right-handed neutrino to decay into a charged
 * lepton, an up-type quark and a down-type quark.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the lepton.
 * @param genq Generation of the quarks.
 * @returns The width and an error estimate.
 */
auto width_vr_to_l_u_d(double mvr, double theta, int genl, int genq) -> double {
  const double ml = kLEPTON_MASSES.at(genl);
  const double mu = kUP_QUARK_MASSES.at(genq);
  const double md = kDOWN_QUARK_MASSES.at(genq);
  // auto msqrd = [mvr, theta, genl,
  //              genq](const std::vector<Pythia8::Vec4> &momenta) {
  //  return msqrd_vr_to_l_u_d(momenta, mvr, theta, genl, genq);
  //};
  // return compute_width_3body(msqrd, mvr, {ml, mu, md}, nevents);

  // auto msqrd = [mvr, theta, genl, genq](double s) {
  //  return partialy_integrated_msqrd_vr_to_l_u_d(s, mvr, theta, genl, genq);
  //};

  PartiallyIntegratedMsqrdParams params{mvr, theta, genl, genq};

  gsl_function integrand;
  integrand.function = &partialy_integrated_msqrd_vr_to_vl_d_d;
  integrand.params = &params;

  return compute_width_3body_quad(&integrand, mvr, {ml, mu, md}).first;
}

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino, and two charged leptons with different generations than neutrino.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the neutrino.
 * @param genlp Generation of charged leptons.
 * @returns The width.
 */
auto width_vr_to_vl_lp_lp(double mvr, double theta, int genl, int genlp)
    -> double {

  const double mlp = kLEPTON_MASSES.at(genlp);

  PartiallyIntegratedMsqrdParams params{mvr, theta, genl, genlp};

  gsl_function integrand;
  integrand.function = &partialy_integrated_msqrd_vr_to_vl_lp_lp;
  integrand.params = &params;

  return compute_width_3body_quad(&integrand, mvr, {0.0, mlp, mlp}).first;
}

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino, and two charged leptons with different generations than neutrino.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the neutrino.
 * @param genlp Generation of charged leptons.
 * @returns The width.
 */
auto width_vr_to_vlp_lp_l(double mvr, double theta, int genl, int genlp)
    -> double {

  const double ml = kLEPTON_MASSES.at(genl);
  const double mlp = kLEPTON_MASSES.at(genlp);

  PartiallyIntegratedMsqrdParams params{mvr, theta, genl, genlp};

  gsl_function integrand;
  integrand.function = &partialy_integrated_msqrd_vr_to_vlp_lp_l;
  integrand.params = &params;

  return compute_width_3body_quad(&integrand, mvr, {0.0, mlp, ml}).first;
}

/**
 * Compute the partial width for a right-handed neutrino to decay into an active
 * neutrino, and two charged leptons with of the same generation as neutrino.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the final-state leptons.
 * @returns The width.
 */
auto width_vr_to_vl_l_l(double mvr, double theta, int genl) -> double {

  const double ml = kLEPTON_MASSES.at(genl);

  PartiallyIntegratedMsqrdParams params{mvr, theta, genl, -1};

  gsl_function integrand;
  integrand.function = &partialy_integrated_msqrd_vr_to_vl_l_l;
  integrand.params = &params;

  return compute_width_3body_quad(&integrand, mvr, {0.0, ml, ml}).first;
}

/**
 * Compute the partial width for a right-handed neutrino to decay into a three
 * active neutrinos, all of the same generation.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the final-state leptons.
 * @returns The width.
 */
auto width_vr_to_vl_vl_vl(double mvr, double theta, int /*genl*/) -> double {

  PartiallyIntegratedMsqrdParams params{mvr, theta, -1, -1};

  gsl_function integrand;
  integrand.function = &partialy_integrated_msqrd_vr_to_vl_vl_vl;
  integrand.params = &params;

  // Extra factor of 1 / 6 us for identical final state particles.
  return compute_width_3body_quad(&integrand, mvr, {0.0, 0.0, 0.0}).first / 6.0;
}

} // namespace storm
