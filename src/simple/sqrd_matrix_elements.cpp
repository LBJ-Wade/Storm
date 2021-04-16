/**
 * This file contains definitions of the squared matrix elements for the decay
 * of a right handed neutrino in the `SimpleRhNeutrino` model.
 */

#include "storm/constants.hpp"
#include "storm/simple.hpp"
#include "storm/types.hpp"
#include <Pythia8/Pythia.h>

namespace storm {

auto SimpleRhNeutrino::msqrd_vr_to_vl_u_u(const MomentaList &momenta) const
    -> double {
  const double s = (momenta.at(1) + momenta.at(2)).m2Calc();
  const double t = (momenta.at(0) + momenta.at(2)).m2Calc();
  const double mu = momenta.at(1).mCalc();

  return (2 * pow(M_PI, 2) * pow(kALPHA_EM, 2) *
          ((-9 * pow(mu, 2) * pow(p_mvr, 2) * (4 * pow(mu, 2) - s) *
            (pow(p_mvr, 2) - s)) /
               (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
           (2 * pow(mu, 4) * pow(kZ_BOSON_MASS, 4) *
                (9 * pow(kCOS_THETA_WEAK, 4) -
                 6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                 17 * pow(kSIN_THETA_WEAK, 4)) +
            pow(kZ_BOSON_MASS, 4) *
                (9 * pow(kCOS_THETA_WEAK, 4) -
                 6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                 17 * pow(kSIN_THETA_WEAK, 4)) *
                (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                 pow(p_mvr, 2) * (s + 2 * t)) +
            pow(mu, 2) *
                (-9 * pow(p_mvr, 4) * (2 * pow(kZ_BOSON_MASS, 2) + s) +
                 9 * pow(p_mvr, 2) *
                     (2 * pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) +
                      2 * pow(kZ_BOSON_MASS, 2) * s + pow(s, 2) +
                      4 * pow(kCOS_THETA_WEAK, 2) * pow(kZ_BOSON_MASS, 4) *
                          pow(kSIN_THETA_WEAK, 2) +
                      2 * pow(kZ_BOSON_MASS, 4) * pow(kSIN_THETA_WEAK, 4)) -
                 2 * pow(kZ_BOSON_MASS, 4) *
                     (6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) *
                          (3 * s - 2 * t) +
                      9 * pow(kCOS_THETA_WEAK, 4) * (s + 2 * t) +
                      pow(kSIN_THETA_WEAK, 4) * (9 * s + 34 * t)))) /
               (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2)))) *
          pow(p_theta, 2)) /
         (pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) *
          pow(kSIN_THETA_WEAK, 4));
}

auto SimpleRhNeutrino::msqrd_vr_to_vl_d_d(const MomentaList &momenta) const
    -> double {

  const double s = (momenta.at(1) + momenta.at(2)).m2Calc();
  const double t = (momenta.at(0) + momenta.at(2)).m2Calc();
  const double md = momenta.at(1).mCalc();

  return (2 * pow(M_PI, 2) * pow(kALPHA_EM, 2) *
          ((-9 * pow(md, 2) * pow(p_mvr, 2) * (4 * pow(md, 2) - s) *
            (pow(p_mvr, 2) - s)) /
               (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
           (2 * pow(md, 4) * pow(kZ_BOSON_MASS, 4) *
                (9 * pow(kCOS_THETA_WEAK, 4) +
                 6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                 5 * pow(kSIN_THETA_WEAK, 4)) +
            pow(kZ_BOSON_MASS, 4) *
                (9 * pow(kCOS_THETA_WEAK, 4) +
                 6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                 5 * pow(kSIN_THETA_WEAK, 4)) *
                (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                 pow(p_mvr, 2) * (s + 2 * t)) -
            pow(md, 2) *
                (9 * pow(p_mvr, 4) * (2 * pow(kZ_BOSON_MASS, 2) + s) -
                 9 * pow(p_mvr, 2) *
                     (2 * pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) +
                      2 * pow(kZ_BOSON_MASS, 2) * s + pow(s, 2) +
                      4 * pow(kCOS_THETA_WEAK, 2) * pow(kZ_BOSON_MASS, 4) *
                          pow(kSIN_THETA_WEAK, 2) +
                      2 * pow(kZ_BOSON_MASS, 4) * pow(kSIN_THETA_WEAK, 4)) +
                 2 * pow(kZ_BOSON_MASS, 4) *
                     (9 * pow(kCOS_THETA_WEAK, 4) * (s + 2 * t) +
                      6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) *
                          (3 * s + 2 * t) +
                      pow(kSIN_THETA_WEAK, 4) * (9 * s + 10 * t)))) /
               (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2)))) *
          pow(p_theta, 2)) /
         (pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) *
          pow(kSIN_THETA_WEAK, 4));
}

auto SimpleRhNeutrino::msqrd_vr_to_l_u_d(const MomentaList &momenta,
                                         size_t gen) const -> double {
  const double s = (momenta.at(1) + momenta.at(2)).m2Calc();
  const double t = (momenta.at(0) + momenta.at(2)).m2Calc();
  const double ml = p_ml;
  const double mu = kUP_QUARK_MASSES.at(gen);
  const double md = kDOWN_QUARK_MASSES.at(gen);
  const double cmk2 = std::abs(kCKM.at(gen).at(gen));

  return (18 * pow(M_PI, 2) *
          ((pow(p_mvr, 2) * (pow(ml, 2) + pow(p_mvr, 2) - s) *
            (-pow(md, 4) - pow(mu, 4) + pow(mu, 2) * s +
             pow(md, 2) * (2 * pow(mu, 2) + s))) /
               (pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4)) +
           (4 * pow(p_mvr, 2) *
            (pow(mu, 2) * (pow(mu, 2) + pow(p_mvr, 2)) -
             pow(md, 2) * (pow(ml, 2) + pow(mu, 2) - t))) /
               (pow(kCOS_THETA_WEAK, 2) * pow(kZ_BOSON_MASS, 2)) +
           4 * (pow(md, 2) + pow(p_mvr, 2) - t) *
               (-pow(ml, 2) - pow(mu, 2) + t)) *
          pow(kALPHA_EM, 2) * pow(p_theta, 2) * cmk2) /
         (pow(kSIN_THETA_WEAK, 4) *
          (pow(kW_BOSON_MASS, 4) + pow(s, 2) +
           pow(kW_BOSON_MASS, 2) * (-2 * s + pow(kW_BOSON_WIDTH, 2))));
}

auto SimpleRhNeutrino::msqrd_vr_to_vl_l_l(const MomentaList &momenta) const
    -> double {

  const double s = (momenta.at(1) + momenta.at(2)).m2Calc();
  const double t = (momenta.at(0) + momenta.at(2)).m2Calc();
  const double ml = p_ml;

  return (2 * pow(M_PI, 2) * pow(kALPHA_EM, 2) *
          ((-4 * (-pow(ml, 2) + s + t) *
            (-pow(ml, 2) - pow(p_mvr, 2) + s + t)) /
               (pow(kW_BOSON_MASS, 4) + pow(t, 2) +
                pow(kW_BOSON_MASS, 2) * (-2 * t + pow(kW_BOSON_WIDTH, 2))) -
           (4 * (pow(ml, 2) - t) * (pow(ml, 2) + pow(p_mvr, 2) - t)) /
               (pow(kW_BOSON_MASS, 4) +
                pow(-2 * pow(ml, 2) - pow(p_mvr, 2) + s + t, 2) +
                pow(kW_BOSON_MASS, 2) *
                    (-4 * pow(ml, 2) - 2 * pow(p_mvr, 2) + 2 * s + 2 * t +
                     pow(kW_BOSON_WIDTH, 2))) +
           (2 *
            (pow(kCOS_THETA_WEAK, 2) * (pow(ml, 2) - s - t) *
                 (pow(ml, 2) + pow(p_mvr, 2) - s - t) -
             pow(kSIN_THETA_WEAK, 2) *
                 (pow(ml, 4) - (pow(p_mvr, 2) - s - t) * (s + t) -
                  pow(ml, 2) * (pow(p_mvr, 2) + 2 * t))) *
            (2 * (pow(kZ_BOSON_MASS, 2) - s) * (pow(kW_BOSON_MASS, 2) - t) +
             2 * kW_BOSON_MASS * kZ_BOSON_MASS * kW_BOSON_WIDTH *
                 kZ_BOSON_WIDTH)) /
               (pow(kCOS_THETA_WEAK, 2) *
                (pow(pow(kW_BOSON_MASS, 2) - t, 2) +
                 pow(kW_BOSON_MASS, 2) * pow(kW_BOSON_WIDTH, 2)) *
                (pow(pow(kZ_BOSON_MASS, 2) - s, 2) +
                 pow(kZ_BOSON_MASS, 2) * pow(kZ_BOSON_WIDTH, 2))) +
           (2 *
            (pow(kCOS_THETA_WEAK, 2) * (pow(ml, 2) - t) *
                 (pow(ml, 2) + pow(p_mvr, 2) - t) -
             pow(kSIN_THETA_WEAK, 2) *
                 (pow(ml, 4) + t * (-pow(p_mvr, 2) + t) -
                  pow(ml, 2) * (pow(p_mvr, 2) - 2 * s + 2 * t))) *
            (2 * (pow(kZ_BOSON_MASS, 2) - s) *
                 (-2 * pow(ml, 2) - pow(p_mvr, 2) + pow(kW_BOSON_MASS, 2) + s +
                  t) +
             2 * kW_BOSON_MASS * kZ_BOSON_MASS * kW_BOSON_WIDTH *
                 kZ_BOSON_WIDTH)) /
               (pow(kCOS_THETA_WEAK, 2) *
                (pow(-2 * pow(ml, 2) - pow(p_mvr, 2) + pow(kW_BOSON_MASS, 2) +
                         s + t,
                     2) +
                 pow(kW_BOSON_MASS, 2) * pow(kW_BOSON_WIDTH, 2)) *
                (pow(pow(kZ_BOSON_MASS, 2) - s, 2) +
                 pow(kZ_BOSON_MASS, 2) * pow(kZ_BOSON_WIDTH, 2))) -
           (pow(kCOS_THETA_WEAK, 4) *
                (2 * pow(ml, 4) + pow(s, 2) +
                 2 * pow(ml, 2) * (pow(p_mvr, 2) - s - 2 * t) + 2 * s * t +
                 2 * pow(t, 2) - pow(p_mvr, 2) * (s + 2 * t)) -
            2 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) *
                (2 * pow(ml, 4) + pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                 2 * pow(ml, 2) * (pow(p_mvr, 2) - s + 2 * t) -
                 pow(p_mvr, 2) * (s + 2 * t)) +
            pow(kSIN_THETA_WEAK, 4) *
                (10 * pow(ml, 4) +
                 2 * pow(ml, 2) * (pow(p_mvr, 2) - s - 10 * t) +
                 5 * (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                      pow(p_mvr, 2) * (s + 2 * t)))) /
               (pow(kCOS_THETA_WEAK, 4) *
                (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                 pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2))))) *
          pow(p_theta, 2)) /
         pow(kSIN_THETA_WEAK, 4);
}

auto SimpleRhNeutrino::msqrd_vr_to_vlp_l_lp(const MomentaList &momenta) const
    -> double {

  const double s = (momenta.at(1) + momenta.at(2)).m2Calc();
  const double t = (momenta.at(0) + momenta.at(2)).m2Calc();
  const double ml1 = p_ml;
  const double ml2 = momenta.at(1).mCalc();

  return (8 * pow(M_PI, 2) * (pow(ml2, 2) + pow(p_mvr, 2) - t) *
          (-pow(ml1, 2) + t) * pow(kALPHA_EM, 2) * pow(p_theta, 2)) /
         (pow(kSIN_THETA_WEAK, 4) *
          (pow(kW_BOSON_MASS, 4) +
           pow(pow(ml1, 2) + pow(ml2, 2) + pow(p_mvr, 2) - s - t, 2) +
           pow(kW_BOSON_MASS, 2) *
               (-2 * (pow(ml1, 2) + pow(ml2, 2) + pow(p_mvr, 2) - s - t) +
                pow(kW_BOSON_WIDTH, 2))));
}

auto SimpleRhNeutrino::msqrd_vr_to_vl_lp_lp(const MomentaList &momenta) const
    -> double {

  const double s = (momenta.at(1) + momenta.at(2)).m2Calc();
  const double t = (momenta.at(0) + momenta.at(2)).m2Calc();
  const double ml = momenta.at(1).mCalc();

  return (-2 * pow(M_PI, 2) *
          (pow(kCOS_THETA_WEAK, 4) *
               (2 * pow(ml, 4) + pow(s, 2) +
                2 * pow(ml, 2) * (pow(p_mvr, 2) - s - 2 * t) + 2 * s * t +
                2 * pow(t, 2) - pow(p_mvr, 2) * (s + 2 * t)) -
           2 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) *
               (2 * pow(ml, 4) + pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                2 * pow(ml, 2) * (pow(p_mvr, 2) - s + 2 * t) -
                pow(p_mvr, 2) * (s + 2 * t)) +
           pow(kSIN_THETA_WEAK, 4) *
               (10 * pow(ml, 4) +
                2 * pow(ml, 2) * (pow(p_mvr, 2) - s - 10 * t) +
                5 * (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                     pow(p_mvr, 2) * (s + 2 * t)))) *
          pow(kALPHA_EM, 2) * pow(p_theta, 2)) /
         (pow(kCOS_THETA_WEAK, 4) * pow(kSIN_THETA_WEAK, 4) *
          (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
           pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2))));
}

} // namespace storm