#include "storm/constants.hpp"
#include "storm/simple.hpp"
#include <Pythia8/Basics.h>
#include <tuple>

namespace storm {

// ===============================
// ---- Two-body final states ----
// ===============================

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and Higgs.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 */
auto msqrd_vr_to_vl_h(const std::vector<Pythia8::Vec4> & /*momenta*/,
                      double mvr, double theta) -> double {
  return (pow(theta, 2) * pow(mvr, 2) * (-pow(kHIGGS_MASS, 2) + pow(mvr, 2))) /
         pow(kHIGGS_VEV, 2);
}

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and Z.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 */
auto msqrd_vr_to_vl_z(const std::vector<Pythia8::Vec4> & /*momenta*/,
                      double mvr, double theta) -> double {
  return ((pow(mvr, 2) - pow(kZ_BOSON_MASS, 2)) *
          (pow(mvr, 2) + 2 * pow(kZ_BOSON_MASS, 2)) * M_PI * kALPHA_EM *
          pow(theta, 2)) /
         (pow(kCOS_THETA_WEAK, 2) * pow(kZ_BOSON_MASS, 2) *
          pow(kSIN_THETA_WEAK, 2));
}

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * a charged lepton and W.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param gen Generation of the final state lepton. Should be 0, 1 or 2.
 */
auto msqrd_vr_to_l_w(const std::vector<Pythia8::Vec4> & /*momenta*/, double mvr,
                     double theta, int gen) -> double {
  const double ml = kLEPTON_MASSES.at(gen);
  return ((pow(ml, 4) + pow(mvr, 4) + pow(mvr, 2) * pow(kW_BOSON_MASS, 2) -
           2 * pow(kW_BOSON_MASS, 4) +
           pow(ml, 2) * (-2 * pow(mvr, 2) + pow(kW_BOSON_MASS, 2))) *
          M_PI * kALPHA_EM * pow(theta, 2)) /
         (pow(kW_BOSON_MASS, 2) * pow(kSIN_THETA_WEAK, 2));
}

// =================================
// ---- Three-body final states ----
// =================================

/**
 * Helper function to compute the Mandelstam variables s, t and u from a set of
 * four-momenta.
 */
static auto momenta_to_stu(const std::vector<Pythia8::Vec4> &momenta)
    -> std::tuple<double, double, double> {
  const double s = (momenta[1] + momenta[2]).m2Calc();
  const double t = (momenta[0] + momenta[2]).m2Calc();
  const double u = (momenta[0] + momenta[1]).m2Calc();
  return std::make_tuple(s, t, u);
}

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and two up-type quarks.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param gen Generation of the down-type quarks.
 */
auto msqrd_vr_to_vl_u_u(const std::vector<Pythia8::Vec4> &momenta, double mvr,
                        double theta, int gen) -> double {
  const double mu = kUP_QUARK_MASSES.at(gen);
  const auto [s, t, u] = momenta_to_stu(momenta);

  return (2 * pow(M_PI, 2) * pow(kALPHA_EM, 2) *
          ((-9 * pow(mu, 2) * pow(mvr, 2) * (4 * pow(mu, 2) - s) *
            (pow(mvr, 2) - s)) /
               (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
           (-2 * pow(mu, 4) *
                (18 * pow(mvr, 2) * pow(kZ_BOSON_MASS, 2) +
                 pow(kZ_BOSON_MASS, 4) *
                     (-9 * pow(kCOS_THETA_WEAK, 4) +
                      6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) -
                      17 * pow(kSIN_THETA_WEAK, 4))) +
            pow(kZ_BOSON_MASS, 4) *
                (9 * pow(kCOS_THETA_WEAK, 4) -
                 6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                 17 * pow(kSIN_THETA_WEAK, 4)) *
                (pow(t, 2) + pow(u, 2) - pow(mvr, 2) * (t + u)) +
            pow(mu, 2) * (-9 * pow(mvr, 4) * s +
                          18 * pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) *
                              (pow(mvr, 2) - t - u) +
                          12 * pow(kCOS_THETA_WEAK, 2) * pow(kZ_BOSON_MASS, 4) *
                              pow(kSIN_THETA_WEAK, 2) *
                              (3 * pow(mvr, 2) - 4 * s + t + u) +
                          2 * pow(kZ_BOSON_MASS, 4) * pow(kSIN_THETA_WEAK, 4) *
                              (8 * s - 17 * (t + u)) +
                          9 * pow(mvr, 2) *
                              (pow(s, 2) + 2 * pow(kZ_BOSON_MASS, 2) *
                                               (pow(kZ_BOSON_MASS, 2) *
                                                    pow(kSIN_THETA_WEAK, 4) +
                                                t + u)))) /
               (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2)))) *
          pow(theta, 2)) /
         (pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) *
          pow(kSIN_THETA_WEAK, 4));
}

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and two down-type quarks.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param gen Generation of the down-type quarks.
 */
auto msqrd_vr_to_vl_d_d(const std::vector<Pythia8::Vec4> &momenta, double mvr,
                        double theta, int gen) -> double {
  const double md = kDOWN_QUARK_MASSES.at(gen);
  const auto [s, t, u] = momenta_to_stu(momenta);

  return (2 * pow(M_PI, 2) * pow(kALPHA_EM, 2) *
          ((-9 * pow(md, 2) * pow(mvr, 2) * (4 * pow(md, 2) - s) *
            (pow(mvr, 2) - s)) /
               (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) +
           (2 * pow(md, 4) *
                (18 * pow(mvr, 2) * pow(kZ_BOSON_MASS, 2) -
                 pow(kZ_BOSON_MASS, 4) *
                     (9 * pow(kCOS_THETA_WEAK, 4) +
                      6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                      5 * pow(kSIN_THETA_WEAK, 4))) +
            pow(kZ_BOSON_MASS, 4) *
                (9 * pow(kCOS_THETA_WEAK, 4) +
                 6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                 5 * pow(kSIN_THETA_WEAK, 4)) *
                (-pow(t, 2) - pow(u, 2) + pow(mvr, 2) * (t + u)) +
            pow(md, 2) * (9 * pow(mvr, 4) * s -
                          18 * pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) *
                              (pow(mvr, 2) - t - u) +
                          12 * pow(kCOS_THETA_WEAK, 2) * pow(kZ_BOSON_MASS, 4) *
                              pow(kSIN_THETA_WEAK, 2) *
                              (-3 * pow(mvr, 2) + 2 * s + t + u) +
                          2 * pow(kZ_BOSON_MASS, 4) * pow(kSIN_THETA_WEAK, 4) *
                              (4 * s + 5 * (t + u)) -
                          9 * pow(mvr, 2) *
                              (pow(s, 2) + 2 * pow(kZ_BOSON_MASS, 2) *
                                               (pow(kZ_BOSON_MASS, 2) *
                                                    pow(kSIN_THETA_WEAK, 4) +
                                                t + u)))) /
               (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2)))) *
          pow(theta, 2)) /
         (pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) *
          pow(kSIN_THETA_WEAK, 4));
}

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * a charged lepton, an up-type quark and an anti-down-type quark.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the charged lepton.
 * @param genq Generation of the quarks.
 */
auto msqrd_vr_to_l_u_d(const std::vector<Pythia8::Vec4> &momenta, double mvr,
                       double theta, int genl, int genq) -> double {

  const double ml = kLEPTON_MASSES.at(genl);
  const double mu = kUP_QUARK_MASSES.at(genq);
  const double md = kDOWN_QUARK_MASSES.at(genq);
  const auto [s, t, u] = momenta_to_stu(momenta);

  return (18 * pow(M_PI, 2) *
          ((pow(mvr, 2) * (pow(ml, 2) + pow(mvr, 2) - s) *
            (-pow(md, 4) - pow(mu, 4) + pow(mu, 2) * s +
             pow(md, 2) * (2 * pow(mu, 2) + s))) /
               (pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4)) -
           (4 * pow(mvr, 2) *
            (pow(mu, 2) * (pow(ml, 2) - t) +
             pow(md, 2) * (pow(ml, 2) + 2 * pow(mu, 2) - u))) /
               (pow(kCOS_THETA_WEAK, 2) * pow(kZ_BOSON_MASS, 2)) +
           4 * (pow(md, 2) + pow(mvr, 2) - u) *
               (-pow(ml, 2) - pow(mu, 2) + u)) *
          pow(kALPHA_EM, 2) * pow(theta, 2)) /
         (pow(kSIN_THETA_WEAK, 4) *
          (pow(kW_BOSON_MASS, 4) + pow(s, 2) +
           pow(kW_BOSON_MASS, 2) * (-2 * s + pow(kW_BOSON_WIDTH, 2))));
}

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
auto msqrd_vr_to_vl_lp_lp(const std::vector<Pythia8::Vec4> &momenta, double mvr,
                          double theta, int /*genl*/, int genlp) -> double {

  const double mlp = kLEPTON_MASSES.at(genlp);
  const auto [s, t, u] = momenta_to_stu(momenta);

  return (-2 * pow(M_PI, 2) *
          (pow(kCOS_THETA_WEAK, 4) *
               (2 * pow(mlp, 4) + pow(s, 2) +
                2 * pow(mlp, 2) * (pow(mvr, 2) - s - 2 * t) + 2 * s * t +
                2 * pow(t, 2) - pow(mvr, 2) * (s + 2 * t)) -
           2 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) *
               (2 * pow(mlp, 4) + pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                2 * pow(mlp, 2) * (pow(mvr, 2) - s + 2 * t) -
                pow(mvr, 2) * (s + 2 * t)) +
           pow(kSIN_THETA_WEAK, 4) *
               (10 * pow(mlp, 4) +
                2 * pow(mlp, 2) * (pow(mvr, 2) - s - 10 * t) +
                5 * (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                     pow(mvr, 2) * (s + 2 * t)))) *
          pow(kALPHA_EM, 2) * pow(theta, 2)) /
         (pow(kCOS_THETA_WEAK, 4) * pow(kSIN_THETA_WEAK, 4) *
          (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
           pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2))));
}

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
auto msqrd_vr_to_vlp_lp_l(const std::vector<Pythia8::Vec4> &momenta, double mvr,
                          double theta, int genl, int genlp) -> double {
  const double ml = kLEPTON_MASSES.at(genl);
  const double mlp = kLEPTON_MASSES.at(genlp);
  const auto [s, t, u] = momenta_to_stu(momenta);
  return (8 * pow(M_PI, 2) * (pow(ml, 2) + pow(mvr, 2) - t) *
          (-pow(mlp, 2) + t) * pow(kALPHA_EM, 2) * pow(theta, 2)) /
         (pow(kSIN_THETA_WEAK, 4) *
          (pow(kW_BOSON_MASS, 4) +
           pow(pow(ml, 2) + pow(mlp, 2) + pow(mvr, 2) - s - t, 2) +
           pow(kW_BOSON_MASS, 2) *
               (-2 * (pow(ml, 2) + pow(mlp, 2) + pow(mvr, 2) - s - t) +
                pow(kW_BOSON_WIDTH, 2))));
}

/**
 * Compute the squared matrix element for a right-handed neutrino to decay into
 * an active neutrino and two charged leptons, all of the same generation.
 * neutrino.
 * @param momenta Four-momenta of the final states.
 * @param mvr Mass of the RH neutrino.
 * @param theta Mixing angle between active and RH neutrino.
 * @param genl Generation of the leptons.
 */
auto msqrd_vr_to_vl_l_l(const std::vector<Pythia8::Vec4> &momenta, double mvr,
                        double theta, int genl) -> double {
  const double ml = kLEPTON_MASSES.at(genl);
  const auto [s, t, u] = momenta_to_stu(momenta);
  return (2 * pow(M_PI, 2) * pow(kALPHA_EM, 2) *
          ((-4 * (-pow(ml, 2) + s + t) * (-pow(ml, 2) - pow(mvr, 2) + s + t)) /
               (pow(kW_BOSON_MASS, 4) + pow(t, 2) +
                pow(kW_BOSON_MASS, 2) * (-2 * t + pow(kW_BOSON_WIDTH, 2))) -
           (4 * (pow(ml, 2) - t) * (pow(ml, 2) + pow(mvr, 2) - t)) /
               (pow(kW_BOSON_MASS, 4) +
                pow(-2 * pow(ml, 2) - pow(mvr, 2) + s + t, 2) +
                pow(kW_BOSON_MASS, 2) *
                    (-4 * pow(ml, 2) - 2 * pow(mvr, 2) + 2 * s + 2 * t +
                     pow(kW_BOSON_WIDTH, 2))) +
           (2 *
            (pow(kCOS_THETA_WEAK, 2) * (pow(ml, 2) - s - t) *
                 (pow(ml, 2) + pow(mvr, 2) - s - t) -
             pow(kSIN_THETA_WEAK, 2) *
                 (pow(ml, 4) - (pow(mvr, 2) - s - t) * (s + t) -
                  pow(ml, 2) * (pow(mvr, 2) + 2 * t))) *
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
                 (pow(ml, 2) + pow(mvr, 2) - t) -
             pow(kSIN_THETA_WEAK, 2) *
                 (pow(ml, 4) + t * (-pow(mvr, 2) + t) -
                  pow(ml, 2) * (pow(mvr, 2) - 2 * s + 2 * t))) *
            (2 * (pow(kZ_BOSON_MASS, 2) - s) *
                 (-2 * pow(ml, 2) - pow(mvr, 2) + pow(kW_BOSON_MASS, 2) + s +
                  t) +
             2 * kW_BOSON_MASS * kZ_BOSON_MASS * kW_BOSON_WIDTH *
                 kZ_BOSON_WIDTH)) /
               (pow(kCOS_THETA_WEAK, 2) *
                (pow(-2 * pow(ml, 2) - pow(mvr, 2) + pow(kW_BOSON_MASS, 2) + s +
                         t,
                     2) +
                 pow(kW_BOSON_MASS, 2) * pow(kW_BOSON_WIDTH, 2)) *
                (pow(pow(kZ_BOSON_MASS, 2) - s, 2) +
                 pow(kZ_BOSON_MASS, 2) * pow(kZ_BOSON_WIDTH, 2))) -
           (pow(kCOS_THETA_WEAK, 4) *
                (2 * pow(ml, 4) + pow(s, 2) +
                 2 * pow(ml, 2) * (pow(mvr, 2) - s - 2 * t) + 2 * s * t +
                 2 * pow(t, 2) - pow(mvr, 2) * (s + 2 * t)) -
            2 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) *
                (2 * pow(ml, 4) + pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                 2 * pow(ml, 2) * (pow(mvr, 2) - s + 2 * t) -
                 pow(mvr, 2) * (s + 2 * t)) +
            pow(kSIN_THETA_WEAK, 4) *
                (10 * pow(ml, 4) + 2 * pow(ml, 2) * (pow(mvr, 2) - s - 10 * t) +
                 5 * (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                      pow(mvr, 2) * (s + 2 * t)))) /
               (pow(kCOS_THETA_WEAK, 4) *
                (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                 pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2))))) *
          pow(theta, 2)) /
         pow(kSIN_THETA_WEAK, 4);
}

} // namespace storm