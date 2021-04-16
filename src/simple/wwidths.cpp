#include "storm/constants.hpp"
#include "storm/rambo.hpp"
#include "storm/simple.hpp"
#include "storm/types.hpp"
#include <array>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <cstddef>
#include <functional>
#include <locale>
#include <stdexcept>
#include <utility>

namespace storm {

/**
 * Compute the partial width for a RH neutrino decaying into a active neutrino
 * and a Higgs.
 */
auto SimpleRhNeutrino::width_vl_h() const -> double {
  return ((-pow(kHIGGS_MASS, 2) + pow(p_mvr, 2)) * pow(p_theta, 2) *
          std::abs(pow(kHIGGS_MASS, 2) - pow(p_mvr, 2))) /
         (16. * p_mvr * M_PI * pow(kHIGGS_VEV, 2));
}

/**
 * Compute the partial width for a RH neutrino decaying into a active neutrino
 * and a Z.
 */
auto SimpleRhNeutrino::width_vl_z() const -> double {

  return (((pow(p_mvr, 4) + pow(p_mvr, 2) * pow(kZ_BOSON_MASS, 2) -
            2 * pow(kZ_BOSON_MASS, 4)) *
           kALPHA_EM * pow(p_theta, 2) *
           std::abs(pow(p_mvr, 2) - pow(kZ_BOSON_MASS, 2))) /
          (16.0 * pow(kCOS_THETA_WEAK, 2) * pow(p_mvr, 3) *
           pow(kZ_BOSON_MASS, 2) * pow(kSIN_THETA_WEAK, 2)));
}

/**
 * Compute the partial width for a RH neutrino decaying into a charged lepton
 * and a W.
 */
auto SimpleRhNeutrino::width_l_w() const -> double {
  return (std::sqrt(
              (p_ml - p_mvr - kW_BOSON_MASS) * (p_ml + p_mvr - kW_BOSON_MASS) *
              (p_ml - p_mvr + kW_BOSON_MASS) * (p_ml + p_mvr + kW_BOSON_MASS)) *
          (pow(pow(p_ml, 2) - pow(p_mvr, 2), 2) +
           (pow(p_ml, 2) + pow(p_mvr, 2)) * pow(kW_BOSON_MASS, 2) -
           2 * pow(kW_BOSON_MASS, 4)) *
          kALPHA_EM * pow(p_theta, 2)) /
         (16. * pow(p_mvr, 3) * pow(kW_BOSON_MASS, 2) *
          pow(kSIN_THETA_WEAK, 2));
}

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
 * @breif Compute the three body decay width given the squared matrix element
 * integrated over the Mandelstam variable `t`.
 * @param msqrd Squared matrix element integrated over t = (p1+p3)^2.
 * @param m Mass of the decaying particle.
 * @param fsp_masses Mass of the three final-state particles.
 * @param error Estimated error in the decay width.
 */
static auto compute_width_3body(const std::function<double(double)> &msqrd,
                                double m, std::array<double, 3> fsp_masses,
                                double *error = nullptr) -> double {
  using boost::math::quadrature::gauss_kronrod;

  auto bounds =
      compute_s_bounds(m, fsp_masses[0], fsp_masses[1], fsp_masses[2]);
  auto integral = gauss_kronrod<double, 15>::integrate(
      msqrd, bounds.first, bounds.second, 5, 1e-9, error);

  const double pf = 1.0 / (256.0 * pow(M_PI * m, 3));
  *error *= pf;
  return std::abs(integral) * pf;
}

static auto compute_width_3body(const SquaredMatrixElement &msqrd,
                                const double m, std::vector<double> fsp_masses,
                                size_t num_events = 10000,
                                double *error = nullptr) {

  Rambo rambo{std::move(fsp_masses), m, msqrd};
  const auto width = rambo.compute_width(num_events);
  *error = width.second;
  return width.first;
}

} // namespace storm