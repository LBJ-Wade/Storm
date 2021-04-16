#include "storm/constants.hpp"
#include "storm/simple.hpp"

namespace storm {

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and up-type quarks.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vl_u_u(double s, void *params) -> double {
  auto *pars = static_cast<PartiallyIntegratedMsqrdParams *>(params);
  const double mvr = pars->mvr;
  const double theta = pars->theta;
  const int genq = pars->genq;
  const double mu = kUP_QUARK_MASSES.at(genq);

  return (-16 * pow(M_PI, 2) *
              (2 * pow(mu, 2) * s + pow(mvr, 2) * s - pow(s, 2) +
               sqrt(s * (-4 * pow(mu, 2) + s) * pow(-pow(mvr, 2) + s, 2))) *
              pow(kALPHA_EM, 2) *
              (3 * pow(kZ_BOSON_MASS, 4) * (2 * pow(mu, 2) + pow(mvr, 2) - s) *
                   s *
                   (2 * pow(mu, 2) * s + pow(mvr, 2) * s - pow(s, 2) +
                    sqrt(s * (-4 * pow(mu, 2) + s) *
                         pow(-pow(mvr, 2) + s, 2))) *
                   (9 * pow(kCOS_THETA_WEAK, 4) -
                    6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                    17 * pow(kSIN_THETA_WEAK, 4)) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               pow(kZ_BOSON_MASS, 4) *
                   pow(2 * pow(mu, 2) * s + pow(mvr, 2) * s - pow(s, 2) +
                           sqrt(s * (-4 * pow(mu, 2) + s) *
                                pow(-pow(mvr, 2) + s, 2)),
                       2) *
                   (9 * pow(kCOS_THETA_WEAK, 4) -
                    6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                    17 * pow(kSIN_THETA_WEAK, 4)) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               6 * pow(s, 2) *
                   (2 * pow(mu, 4) * pow(kZ_BOSON_MASS, 4) *
                        (9 * pow(kCOS_THETA_WEAK, 4) -
                         6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                         17 * pow(kSIN_THETA_WEAK, 4)) -
                    pow(kZ_BOSON_MASS, 4) * (pow(mvr, 2) - s) * s *
                        (9 * pow(kCOS_THETA_WEAK, 4) -
                         6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                         17 * pow(kSIN_THETA_WEAK, 4)) +
                    9 * pow(mu, 2) *
                        (pow(mvr, 4) * (2 * pow(kZ_BOSON_MASS, 2) - s) -
                         2 * pow(kZ_BOSON_MASS, 4) * s *
                             pow(pow(kCOS_THETA_WEAK, 2) +
                                     pow(kSIN_THETA_WEAK, 2),
                                 2) +
                         pow(mvr, 2) *
                             (2 * pow(kCOS_THETA_WEAK, 4) *
                                  pow(kZ_BOSON_MASS, 4) -
                              2 * pow(kZ_BOSON_MASS, 2) * s + pow(s, 2) +
                              4 * pow(kCOS_THETA_WEAK, 2) *
                                  pow(kZ_BOSON_MASS, 4) *
                                  pow(kSIN_THETA_WEAK, 2) +
                              2 * pow(kZ_BOSON_MASS, 4) *
                                  pow(kSIN_THETA_WEAK, 4)))) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               54 * pow(mu, 2) * pow(mvr, 2) * (4 * pow(mu, 2) - s) *
                   (pow(mvr, 2) - s) * pow(s, 2) *
                   (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                    pow(kZ_BOSON_MASS, 2) *
                        (-2 * s + pow(kZ_BOSON_WIDTH, 2)))) *
              pow(theta, 2) +
          16 * pow(M_PI, 2) *
              (2 * pow(mu, 2) * s + pow(mvr, 2) * s - pow(s, 2) -
               sqrt(s * (-4 * pow(mu, 2) + s) * pow(-pow(mvr, 2) + s, 2))) *
              pow(kALPHA_EM, 2) *
              (3 * pow(kZ_BOSON_MASS, 4) * (2 * pow(mu, 2) + pow(mvr, 2) - s) *
                   s *
                   (2 * pow(mu, 2) * s + pow(mvr, 2) * s - pow(s, 2) -
                    sqrt(s * (-4 * pow(mu, 2) + s) *
                         pow(-pow(mvr, 2) + s, 2))) *
                   (9 * pow(kCOS_THETA_WEAK, 4) -
                    6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                    17 * pow(kSIN_THETA_WEAK, 4)) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               pow(kZ_BOSON_MASS, 4) *
                   pow(-2 * pow(mu, 2) * s - pow(mvr, 2) * s + pow(s, 2) +
                           sqrt(s * (-4 * pow(mu, 2) + s) *
                                pow(-pow(mvr, 2) + s, 2)),
                       2) *
                   (9 * pow(kCOS_THETA_WEAK, 4) -
                    6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                    17 * pow(kSIN_THETA_WEAK, 4)) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               6 * pow(s, 2) *
                   (2 * pow(mu, 4) * pow(kZ_BOSON_MASS, 4) *
                        (9 * pow(kCOS_THETA_WEAK, 4) -
                         6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                         17 * pow(kSIN_THETA_WEAK, 4)) -
                    pow(kZ_BOSON_MASS, 4) * (pow(mvr, 2) - s) * s *
                        (9 * pow(kCOS_THETA_WEAK, 4) -
                         6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                         17 * pow(kSIN_THETA_WEAK, 4)) +
                    9 * pow(mu, 2) *
                        (pow(mvr, 4) * (2 * pow(kZ_BOSON_MASS, 2) - s) -
                         2 * pow(kZ_BOSON_MASS, 4) * s *
                             pow(pow(kCOS_THETA_WEAK, 2) +
                                     pow(kSIN_THETA_WEAK, 2),
                                 2) +
                         pow(mvr, 2) *
                             (2 * pow(kCOS_THETA_WEAK, 4) *
                                  pow(kZ_BOSON_MASS, 4) -
                              2 * pow(kZ_BOSON_MASS, 2) * s + pow(s, 2) +
                              4 * pow(kCOS_THETA_WEAK, 2) *
                                  pow(kZ_BOSON_MASS, 4) *
                                  pow(kSIN_THETA_WEAK, 2) +
                              2 * pow(kZ_BOSON_MASS, 4) *
                                  pow(kSIN_THETA_WEAK, 4)))) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               54 * pow(mu, 2) * pow(mvr, 2) * (4 * pow(mu, 2) - s) *
                   (pow(mvr, 2) - s) * pow(s, 2) *
                   (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                    pow(kZ_BOSON_MASS, 2) *
                        (-2 * s + pow(kZ_BOSON_WIDTH, 2)))) *
              pow(theta, 2)) /
         (96. * pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) * pow(s, 3) *
          pow(kSIN_THETA_WEAK, 4) *
          (pow(kHIGGS_MASS, 4) + pow(s, 2) +
           pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) *
          (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
           pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2))));
}

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and down-type quarks.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vl_d_d(double s, void *params) -> double {
  auto *pars = static_cast<PartiallyIntegratedMsqrdParams *>(params);
  const double mvr = pars->mvr;
  const double theta = pars->theta;
  const int genq = pars->genq;
  const double md = kDOWN_QUARK_MASSES.at(genq);

  return (-16 * pow(M_PI, 2) *
              (2 * pow(md, 2) * s + pow(mvr, 2) * s - pow(s, 2) +
               sqrt(s * (-4 * pow(md, 2) + s) * pow(-pow(mvr, 2) + s, 2))) *
              pow(kALPHA_EM, 2) *
              (3 * pow(kZ_BOSON_MASS, 4) * (2 * pow(md, 2) + pow(mvr, 2) - s) *
                   s *
                   (2 * pow(md, 2) * s + pow(mvr, 2) * s - pow(s, 2) +
                    sqrt(s * (-4 * pow(md, 2) + s) *
                         pow(-pow(mvr, 2) + s, 2))) *
                   (9 * pow(kCOS_THETA_WEAK, 4) +
                    6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                    5 * pow(kSIN_THETA_WEAK, 4)) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               pow(kZ_BOSON_MASS, 4) *
                   pow(2 * pow(md, 2) * s + pow(mvr, 2) * s - pow(s, 2) +
                           sqrt(s * (-4 * pow(md, 2) + s) *
                                pow(-pow(mvr, 2) + s, 2)),
                       2) *
                   (9 * pow(kCOS_THETA_WEAK, 4) +
                    6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                    5 * pow(kSIN_THETA_WEAK, 4)) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               6 * pow(s, 2) *
                   (2 * pow(md, 4) * pow(kZ_BOSON_MASS, 4) *
                        (9 * pow(kCOS_THETA_WEAK, 4) +
                         6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                         5 * pow(kSIN_THETA_WEAK, 4)) -
                    pow(kZ_BOSON_MASS, 4) * (pow(mvr, 2) - s) * s *
                        (9 * pow(kCOS_THETA_WEAK, 4) +
                         6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                         5 * pow(kSIN_THETA_WEAK, 4)) +
                    9 * pow(md, 2) *
                        (pow(mvr, 4) * (2 * pow(kZ_BOSON_MASS, 2) - s) -
                         2 * pow(kZ_BOSON_MASS, 4) * s *
                             pow(pow(kCOS_THETA_WEAK, 2) +
                                     pow(kSIN_THETA_WEAK, 2),
                                 2) +
                         pow(mvr, 2) *
                             (2 * pow(kCOS_THETA_WEAK, 4) *
                                  pow(kZ_BOSON_MASS, 4) -
                              2 * pow(kZ_BOSON_MASS, 2) * s + pow(s, 2) +
                              4 * pow(kCOS_THETA_WEAK, 2) *
                                  pow(kZ_BOSON_MASS, 4) *
                                  pow(kSIN_THETA_WEAK, 2) +
                              2 * pow(kZ_BOSON_MASS, 4) *
                                  pow(kSIN_THETA_WEAK, 4)))) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               54 * pow(md, 2) * pow(mvr, 2) * (4 * pow(md, 2) - s) *
                   (pow(mvr, 2) - s) * pow(s, 2) *
                   (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                    pow(kZ_BOSON_MASS, 2) *
                        (-2 * s + pow(kZ_BOSON_WIDTH, 2)))) *
              pow(theta, 2) +
          16 * pow(M_PI, 2) *
              (2 * pow(md, 2) * s + pow(mvr, 2) * s - pow(s, 2) -
               sqrt(s * (-4 * pow(md, 2) + s) * pow(-pow(mvr, 2) + s, 2))) *
              pow(kALPHA_EM, 2) *
              (3 * pow(kZ_BOSON_MASS, 4) * (2 * pow(md, 2) + pow(mvr, 2) - s) *
                   s *
                   (2 * pow(md, 2) * s + pow(mvr, 2) * s - pow(s, 2) -
                    sqrt(s * (-4 * pow(md, 2) + s) *
                         pow(-pow(mvr, 2) + s, 2))) *
                   (9 * pow(kCOS_THETA_WEAK, 4) +
                    6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                    5 * pow(kSIN_THETA_WEAK, 4)) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               pow(kZ_BOSON_MASS, 4) *
                   pow(-2 * pow(md, 2) * s - pow(mvr, 2) * s + pow(s, 2) +
                           sqrt(s * (-4 * pow(md, 2) + s) *
                                pow(-pow(mvr, 2) + s, 2)),
                       2) *
                   (9 * pow(kCOS_THETA_WEAK, 4) +
                    6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                    5 * pow(kSIN_THETA_WEAK, 4)) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               6 * pow(s, 2) *
                   (2 * pow(md, 4) * pow(kZ_BOSON_MASS, 4) *
                        (9 * pow(kCOS_THETA_WEAK, 4) +
                         6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                         5 * pow(kSIN_THETA_WEAK, 4)) -
                    pow(kZ_BOSON_MASS, 4) * (pow(mvr, 2) - s) * s *
                        (9 * pow(kCOS_THETA_WEAK, 4) +
                         6 * pow(kCOS_THETA_WEAK, 2) * pow(kSIN_THETA_WEAK, 2) +
                         5 * pow(kSIN_THETA_WEAK, 4)) +
                    9 * pow(md, 2) *
                        (pow(mvr, 4) * (2 * pow(kZ_BOSON_MASS, 2) - s) -
                         2 * pow(kZ_BOSON_MASS, 4) * s *
                             pow(pow(kCOS_THETA_WEAK, 2) +
                                     pow(kSIN_THETA_WEAK, 2),
                                 2) +
                         pow(mvr, 2) *
                             (2 * pow(kCOS_THETA_WEAK, 4) *
                                  pow(kZ_BOSON_MASS, 4) -
                              2 * pow(kZ_BOSON_MASS, 2) * s + pow(s, 2) +
                              4 * pow(kCOS_THETA_WEAK, 2) *
                                  pow(kZ_BOSON_MASS, 4) *
                                  pow(kSIN_THETA_WEAK, 2) +
                              2 * pow(kZ_BOSON_MASS, 4) *
                                  pow(kSIN_THETA_WEAK, 4)))) *
                   (pow(kHIGGS_MASS, 4) + pow(s, 2) +
                    pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) -
               54 * pow(md, 2) * pow(mvr, 2) * (4 * pow(md, 2) - s) *
                   (pow(mvr, 2) - s) * pow(s, 2) *
                   (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
                    pow(kZ_BOSON_MASS, 2) *
                        (-2 * s + pow(kZ_BOSON_WIDTH, 2)))) *
              pow(theta, 2)) /
         (96. * pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) * pow(s, 3) *
          pow(kSIN_THETA_WEAK, 4) *
          (pow(kHIGGS_MASS, 4) + pow(s, 2) +
           pow(kHIGGS_MASS, 2) * (-2 * s + pow(kHIGGS_WIDTH, 2))) *
          (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
           pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2))));
}

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into a charged lepton, an up-type quark and a down-type
 * quark.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_l_u_d(double s, void *params) -> double {

  auto *pars = static_cast<PartiallyIntegratedMsqrdParams *>(params);
  const double mvr = pars->mvr;
  const double theta = pars->theta;
  const int genl = pars->genl;
  const int genq = pars->genq;

  const double ml = kLEPTON_MASSES.at(genl);
  const double mu = kUP_QUARK_MASSES.at(genq);
  const double md = kDOWN_QUARK_MASSES.at(genq);

  return (-6 * pow(M_PI, 2) *
          sqrt(pow(md, 4) + pow(pow(mu, 2) - s, 2) -
               2 * pow(md, 2) * (pow(mu, 2) + s)) *
          sqrt(pow(ml, 4) + pow(pow(mvr, 2) - s, 2) -
               2 * pow(ml, 2) * (pow(mvr, 2) + s)) *
          (-3 * pow(mvr, 2) * (pow(ml, 2) + pow(mvr, 2) - s) * pow(s, 2) *
               (pow(md, 4) + pow(mu, 4) - pow(mu, 2) * s -
                pow(md, 2) * (2 * pow(mu, 2) + s)) +
           6 * pow(kCOS_THETA_WEAK, 2) * pow(mvr, 2) * pow(kZ_BOSON_MASS, 2) *
               s * (pow(ml, 2) - pow(mvr, 2) + s) *
               (pow(md, 4) + pow(mu, 4) - pow(mu, 2) * s -
                pow(md, 2) * (2 * pow(mu, 2) + s)) +
           2 * pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) *
               (pow(md, 4) *
                    (-2 * pow(ml, 4) - 2 * pow(mvr, 4) + pow(mvr, 2) * s +
                     pow(s, 2) + pow(ml, 2) * (4 * pow(mvr, 2) + s)) +
                pow(md, 2) * (pow(ml, 4) * (4 * pow(mu, 2) + s) +
                              (pow(mvr, 2) - s) *
                                  ((pow(mvr, 2) - s) * s +
                                   2 * pow(mu, 2) * (2 * pow(mvr, 2) + s)) -
                              2 * pow(ml, 2) *
                                  (s * (pow(mvr, 2) + s) +
                                   pow(mu, 2) * (4 * pow(mvr, 2) + s))) +
                (pow(mu, 2) - s) *
                    (-(pow(ml, 4) * (2 * pow(mu, 2) + s)) +
                     pow(ml, 2) * ((2 * pow(mvr, 2) - s) * s +
                                   pow(mu, 2) * (4 * pow(mvr, 2) + s)) -
                     (pow(mvr, 2) - s) * (pow(mu, 2) * (2 * pow(mvr, 2) + s) +
                                          s * (pow(mvr, 2) + 2 * s))))) *
          pow(kALPHA_EM, 2) * pow(theta, 2)) /
         (pow(kCOS_THETA_WEAK, 4) * pow(kZ_BOSON_MASS, 4) * pow(s, 3) *
          pow(kSIN_THETA_WEAK, 4) *
          (pow(kW_BOSON_MASS, 4) + pow(s, 2) +
           pow(kW_BOSON_MASS, 2) * (-2 * s + pow(kW_BOSON_WIDTH, 2))));
}

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and two charged leptons of
 * different generation than neutrinos.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vl_lp_lp(double s, void *params)
    -> double {
  auto *pars = static_cast<PartiallyIntegratedMsqrdParams *>(params);
  const double mvr = pars->mvr;
  const double theta = pars->theta;
  const int genlp = pars->genq;

  const double mlp = kLEPTON_MASSES.at(genlp);
  return (-2 * pow(M_PI, 2) * (pow(mvr, 2) - s) *
          sqrt(pow(pow(mvr, 2) - s, 2) * s * (-4 * pow(mlp, 2) + s)) *
          (pow(kCOS_THETA_WEAK, 4) * (2 * pow(mlp, 2) * (pow(mvr, 2) - s) +
                                      s * (pow(mvr, 2) + 2 * s)) -
           2 * pow(kCOS_THETA_WEAK, 2) *
               (s * (pow(mvr, 2) + 2 * s) +
                2 * pow(mlp, 2) * (pow(mvr, 2) + 5 * s)) *
               pow(kSIN_THETA_WEAK, 2) +
           (5 * s * (pow(mvr, 2) + 2 * s) +
            2 * pow(mlp, 2) * (5 * pow(mvr, 2) + 7 * s)) *
               pow(kSIN_THETA_WEAK, 4)) *
          pow(kALPHA_EM, 2) * pow(theta, 2)) /
         (3. * pow(kCOS_THETA_WEAK, 4) * pow(s, 2) * pow(kSIN_THETA_WEAK, 4) *
          (pow(kZ_BOSON_MASS, 4) + pow(s, 2) +
           pow(kZ_BOSON_MASS, 2) * (-2 * s + pow(kZ_BOSON_WIDTH, 2))));
}

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and two charged leptons. Here
 * the neutrino shares generation with the charged lepton and the right-haned
 * neutrino is of the same generation as anti-charged lepton.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vlp_lp_l(double s, void *params)
    -> double {

  auto *pars = static_cast<PartiallyIntegratedMsqrdParams *>(params);
  const double mvr = pars->mvr;
  const double theta = pars->theta;
  const int genl = pars->genl;
  const int genq = pars->genq;

  const double ml = kLEPTON_MASSES.at(genl);
  const double mlp = kLEPTON_MASSES.at(genq);

  return (4 * pow(M_PI, 2) * pow(kALPHA_EM, 2) * pow(theta, 2) *
          ((pow(ml, 2) * (pow(mvr, 2) + s) +
            (-pow(mvr, 2) + s) *
                (pow(mlp, 2) - s +
                 std::sqrt(pow(ml, 4) + pow(pow(mlp, 2) - s, 2) -
                           2 * pow(ml, 2) * (pow(mlp, 2) + s)))) /
               s +
           (-(pow(ml, 2) * (pow(mvr, 2) + s)) +
            (-pow(mvr, 2) + s) *
                (-pow(mlp, 2) + s +
                 std::sqrt(pow(ml, 4) + pow(pow(mlp, 2) - s, 2) -
                           2 * pow(ml, 2) * (pow(mlp, 2) + s)))) /
               s -
           (2 *
            (-(pow(mvr, 2) * pow(kW_BOSON_MASS, 2)) + pow(kW_BOSON_MASS, 4) +
             pow(ml, 2) * (pow(mlp, 2) - pow(kW_BOSON_MASS, 2) - s) +
             pow(mlp, 2) * (pow(mvr, 2) - pow(kW_BOSON_MASS, 2) - s) -
             pow(mvr, 2) * s + 2 * pow(kW_BOSON_MASS, 2) * s + pow(s, 2) -
             pow(kW_BOSON_MASS, 2) * pow(kW_BOSON_WIDTH, 2)) *
            std::atan((pow(ml, 2) + pow(mlp, 2) + pow(mvr, 2) -
                       pow(kW_BOSON_MASS, 2) - s -
                       (pow(ml, 2) * (pow(mvr, 2) + s) +
                        (-pow(mvr, 2) + s) *
                            (pow(mlp, 2) - s +
                             std::sqrt(pow(ml, 4) + pow(pow(mlp, 2) - s, 2) -
                                       2 * pow(ml, 2) * (pow(mlp, 2) + s)))) /
                           (2. * s)) /
                      (kW_BOSON_MASS * kW_BOSON_WIDTH))) /
               (kW_BOSON_MASS * kW_BOSON_WIDTH) +
           (2 *
            (-(pow(mvr, 2) * pow(kW_BOSON_MASS, 2)) + pow(kW_BOSON_MASS, 4) +
             pow(ml, 2) * (pow(mlp, 2) - pow(kW_BOSON_MASS, 2) - s) +
             pow(mlp, 2) * (pow(mvr, 2) - pow(kW_BOSON_MASS, 2) - s) -
             pow(mvr, 2) * s + 2 * pow(kW_BOSON_MASS, 2) * s + pow(s, 2) -
             pow(kW_BOSON_MASS, 2) * pow(kW_BOSON_WIDTH, 2)) *
            std::atan((pow(ml, 2) + pow(mlp, 2) + pow(mvr, 2) -
                       pow(kW_BOSON_MASS, 2) - s -
                       (pow(ml, 2) * (pow(mvr, 2) + s) +
                        (pow(mvr, 2) - s) *
                            (-pow(mlp, 2) + s +
                             std::sqrt(pow(ml, 4) + pow(pow(mlp, 2) - s, 2) -
                                       2 * pow(ml, 2) * (pow(mlp, 2) + s)))) /
                           (2. * s)) /
                      (kW_BOSON_MASS * kW_BOSON_WIDTH))) /
               (kW_BOSON_MASS * kW_BOSON_WIDTH) +
           (pow(ml, 2) + pow(mlp, 2) + pow(mvr, 2) - 2 * pow(kW_BOSON_MASS, 2) -
            2 * s) *
               std::log(
                   (pow(kW_BOSON_MASS, 4) +
                    (pow(kW_BOSON_MASS, 2) *
                     (pow(ml, 2) * (pow(mvr, 2) - s) -
                      pow(mlp, 2) * (pow(mvr, 2) + s) +
                      (-pow(mvr, 2) + s) *
                          (s +
                           std::sqrt(pow(ml, 4) + pow(pow(mlp, 2) - s, 2) -
                                     2 * pow(ml, 2) * (pow(mlp, 2) + s))))) /
                        s +
                    pow(pow(ml, 2) * (pow(mvr, 2) - s) -
                            pow(mlp, 2) * (pow(mvr, 2) + s) +
                            (-pow(mvr, 2) + s) *
                                (s + std::sqrt(
                                         pow(ml, 4) + pow(pow(mlp, 2) - s, 2) -
                                         2 * pow(ml, 2) * (pow(mlp, 2) + s))),
                        2) /
                        (4. * pow(s, 2)) +
                    pow(kW_BOSON_MASS, 2) * pow(kW_BOSON_WIDTH, 2)) /
                   (pow(kW_BOSON_MASS, 4) -
                    (pow(kW_BOSON_MASS, 2) *
                     (pow(ml, 2) * (-pow(mvr, 2) + s) +
                      pow(mlp, 2) * (pow(mvr, 2) + s) +
                      (-pow(mvr, 2) + s) *
                          (-s +
                           std::sqrt(pow(ml, 4) + pow(pow(mlp, 2) - s, 2) -
                                     2 * pow(ml, 2) * (pow(mlp, 2) + s))))) /
                        s +
                    pow(pow(ml, 2) * (-pow(mvr, 2) + s) +
                            pow(mlp, 2) * (pow(mvr, 2) + s) +
                            (-pow(mvr, 2) + s) *
                                (-s + std::sqrt(
                                          pow(ml, 4) + pow(pow(mlp, 2) - s, 2) -
                                          2 * pow(ml, 2) * (pow(mlp, 2) + s))),
                        2) /
                        (4. * pow(s, 2)) +
                    pow(kW_BOSON_MASS, 2) * pow(kW_BOSON_WIDTH, 2))))) /
         pow(kSIN_THETA_WEAK, 4);
}

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an active neutrino and two charged lepton, all
 * with the same generation. quark.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta, genl and genq.
 */
auto partialy_integrated_msqrd_vr_to_vl_l_l(double s, void *params) -> double {

  auto *pars = static_cast<PartiallyIntegratedMsqrdParams *>(params);
  const double mvr = pars->mvr;
  const double theta = pars->theta;
  const int genl = pars->genl;
  const double ml = kLEPTON_MASSES.at(genl);

  constexpr double temp1 = kZ_BOSON_MASS * kZ_BOSON_MASS;
  double temp2 = pow(s, 2);
  double temp3 = pow(ml, 2);
  double temp4 = pow(mvr, 2);
  double temp5 = pow(kCOS_THETA_WEAK, 2);
  double temp6 = -s;
  double temp7 = temp1 + temp6;
  double temp8 = pow(s, -2);
  double temp9 = pow(kW_BOSON_MASS, 2);
  double temp10 = -4 * temp3;
  double temp11 = s + temp10;
  double temp12 = -temp4;
  double temp13 = s + temp12;
  double temp14 = pow(temp13, 2);
  double temp15 = s * temp11 * temp14;
  double temp16 = sqrt(temp15);
  double temp17 = pow(kSIN_THETA_WEAK, 2);
  double temp18 = -temp17;
  double temp19 = temp18 + temp5;
  double temp20 = 2 * s * temp3;
  double temp21 = s * temp4;
  double temp22 = -temp2;
  double temp23 = temp16 + temp20 + temp21 + temp22;
  double temp24 = pow(kCOS_THETA_WEAK, 4);
  double temp25 = -2 * temp17 * temp5;
  double temp26 = pow(kSIN_THETA_WEAK, 4);
  double temp27 = 5 * temp26;
  double temp28 = temp24 + temp25 + temp27;
  double temp29 = 2 * temp3;
  double temp30 = temp29 + temp4 + temp6;
  double temp31 = -2 * s * temp3;
  double temp32 = -(s * temp4);
  double temp33 = pow(s, -3);
  double temp34 = temp16 + temp2 + temp31 + temp32;
  double temp35 = pow(ml, 4);
  double temp36 = 2 * temp35;
  double temp37 = temp4 + temp6;
  double temp38 = s * temp13;
  double temp39 = 2 * temp3 * temp37;
  double temp40 = 1 / s;
  double temp41 = temp36 + temp38 + temp39;
  double temp42 = temp24 * temp41;
  double temp43 = -2 * temp3 * temp37;
  double temp44 = temp36 + temp38 + temp43;
  double temp45 = -2 * temp17 * temp44 * temp5;
  double temp46 = 10 * temp35;
  double temp47 = 5 * s * temp13;
  double temp48 = temp39 + temp46 + temp47;
  double temp49 = temp26 * temp48;
  double temp50 = temp42 + temp45 + temp49;
  double temp51 = 2 * s * temp9;
  double temp52 = temp16 + temp2 + temp31 + temp32 + temp51;
  double temp53 = -2 * temp1 * temp9;
  double temp54 = 2 * temp3 * temp7;
  double temp55 = temp4 * temp7;
  double temp56 = -2 * s * temp1;
  double temp57 = 2 * temp2;
  double temp58 =
      kW_BOSON_MASS * kZ_BOSON_MASS * kW_BOSON_WIDTH * kZ_BOSON_WIDTH;
  double temp59 = temp51 + temp53 + temp54 + temp55 + temp56 + temp57 + temp58;
  double temp60 = pow(kZ_BOSON_MASS, 4);
  double temp61 = -2 * s;
  double temp62 = pow(kZ_BOSON_WIDTH, 2);
  double temp63 = temp61 + temp62;
  double temp64 = temp1 * temp63;
  double temp65 = temp2 + temp60 + temp64;
  double temp66 = s + temp9;
  double temp67 = 1 / kW_BOSON_MASS;
  double temp68 = -2 * s * temp9;
  double temp69 = temp16 + temp20 + temp21 + temp22 + temp68;
  double temp70 = 1 / kW_BOSON_WIDTH;
  double temp71 = pow(kW_BOSON_MASS, 3);
  double temp72 = pow(kW_BOSON_MASS, 4);
  double temp73 = pow(kW_BOSON_WIDTH, 2);
  double temp74 = kW_BOSON_MASS * temp7 * kW_BOSON_WIDTH;
  double temp75 = kZ_BOSON_MASS * temp9 * kZ_BOSON_WIDTH;
  double temp76 = kZ_BOSON_MASS * s * kZ_BOSON_WIDTH;
  double temp77 = temp74 + temp75 + temp76;
  double temp78 = (temp40 * temp67 * temp69 * temp70) / 2.;
  double temp79 = std::atan(temp78);
  double temp80 = -(temp4 * temp66);
  double temp81 = -2 * temp66;
  double temp82 = temp4 + temp81;
  double temp83 = temp3 * temp82;
  double temp84 = -(temp73 * temp9);
  double temp85 = temp2 + temp35 + temp51 + temp72 + temp80 + temp83 + temp84;
  double temp86 = 2 * temp1 * temp71 * kW_BOSON_WIDTH;
  double temp87 = -2 * s * temp71 * kW_BOSON_WIDTH;
  double temp88 = 2 * kW_BOSON_MASS * s * temp1 * kW_BOSON_WIDTH;
  double temp89 = -2 * kW_BOSON_MASS * temp2 * kW_BOSON_WIDTH;
  double temp90 = kZ_BOSON_MASS * temp35 * kZ_BOSON_WIDTH;
  double temp91 = kZ_BOSON_MASS * temp72 * kZ_BOSON_WIDTH;
  double temp92 = 2 * kZ_BOSON_MASS * s * temp9 * kZ_BOSON_WIDTH;
  double temp93 = kZ_BOSON_MASS * temp2 * kZ_BOSON_WIDTH;
  double temp94 = -(kZ_BOSON_MASS * temp73 * temp9 * kZ_BOSON_WIDTH);
  double temp95 = -temp1;
  double temp96 = s + temp95;
  double temp97 = 2 * kW_BOSON_MASS * temp96 * kW_BOSON_WIDTH;
  double temp98 = -2 * kZ_BOSON_MASS * temp9 * kZ_BOSON_WIDTH;
  double temp99 = temp4 + temp61;
  double temp100 = kZ_BOSON_MASS * temp99 * kZ_BOSON_WIDTH;
  double temp101 = temp100 + temp97 + temp98;
  double temp102 = temp101 * temp3;
  double temp103 = -(temp4 * temp77);
  double temp104 = temp102 + temp103 + temp86 + temp87 + temp88 + temp89 +
                   temp90 + temp91 + temp92 + temp93 + temp94;
  double temp105 = temp104 * temp5;
  double temp106 = -2 * temp1 * temp71 * kW_BOSON_WIDTH;
  double temp107 = 2 * s * temp71 * kW_BOSON_WIDTH;
  double temp108 = -2 * kW_BOSON_MASS * s * temp1 * kW_BOSON_WIDTH;
  double temp109 = 2 * kW_BOSON_MASS * temp2 * kW_BOSON_WIDTH;
  double temp110 = -(kZ_BOSON_MASS * temp35 * kZ_BOSON_WIDTH);
  double temp111 = -(kZ_BOSON_MASS * temp72 * kZ_BOSON_WIDTH);
  double temp112 = -2 * kZ_BOSON_MASS * s * temp9 * kZ_BOSON_WIDTH;
  double temp113 = -(kZ_BOSON_MASS * temp2 * kZ_BOSON_WIDTH);
  double temp114 = kZ_BOSON_MASS * temp73 * temp9 * kZ_BOSON_WIDTH;
  double temp115 = 2 * kW_BOSON_MASS * temp7 * kW_BOSON_WIDTH;
  double temp116 = kZ_BOSON_MASS * temp4 * kZ_BOSON_WIDTH;
  double temp117 = 2 * kZ_BOSON_MASS * temp9 * kZ_BOSON_WIDTH;
  double temp118 = temp115 + temp116 + temp117;
  double temp119 = temp118 * temp3;
  double temp120 = temp4 * temp77;
  double temp121 = temp106 + temp107 + temp108 + temp109 + temp110 + temp111 +
                   temp112 + temp113 + temp114 + temp119 + temp120;
  double temp122 = temp121 * temp17;
  double temp123 = temp105 + temp122;
  double temp124 = (temp40 * temp52 * temp67 * temp70) / 2.;
  double temp125 = std::atan(temp124);
  double temp126 = pow(temp23, 2);
  double temp127 = pow(temp34, 2);
  double temp128 = temp73 * temp9;
  double temp129 = pow(s, 3);
  double temp130 = temp9 * temp96;
  double temp131 = s * temp96;
  double temp132 = temp130 + temp131 + temp58;
  double temp133 = pow(temp69, 2);
  double temp134 = pow(temp52, 2);
  double temp135 = -(temp23 * temp40) / 2.;
  double temp136 = temp135 + temp9;
  double temp137 = (temp34 * temp40) / 2.;
  double temp138 = temp137 + temp9;
  double temp139 = temp29 + temp4 + temp81;
  double temp140 = (temp126 * temp8) / 4.;
  double temp141 = (temp127 * temp8) / 4.;
  return (pow(kQE, 4) * temp67 * temp70 * pow(theta, 2) *
          (24 * temp125 * temp24 * temp65 * temp85 +
           24 * temp24 * temp65 * temp79 * temp85 +
           (kW_BOSON_MASS * pow(temp23, 3) * temp28 * temp33 * kW_BOSON_WIDTH) /
               4. +
           (kW_BOSON_MASS * temp28 * temp33 * pow(temp34, 3) * kW_BOSON_WIDTH) /
               4. -
           12 * kW_BOSON_MASS * temp123 * temp125 * temp5 * kW_BOSON_WIDTH +
           (3 * kW_BOSON_MASS * temp23 * temp40 * temp50 * kW_BOSON_WIDTH) /
               2. +
           (3 * kW_BOSON_MASS * temp34 * temp40 * temp50 * kW_BOSON_WIDTH) /
               2. +
           12 * kW_BOSON_MASS * temp19 * temp40 * temp5 *
               (-temp16 + temp2 + temp31 + temp32 + temp51) * temp59 *
               kW_BOSON_WIDTH -
           12 * kW_BOSON_MASS * temp19 * temp40 * temp5 * temp52 * temp59 *
               kW_BOSON_WIDTH +
           12 * kW_BOSON_MASS * temp23 * temp24 * temp40 * temp65 *
               kW_BOSON_WIDTH +
           12 * kW_BOSON_MASS * temp24 * temp34 * temp40 * temp65 *
               kW_BOSON_WIDTH -
           12 * kW_BOSON_MASS * temp123 * temp5 * temp79 * kW_BOSON_WIDTH -
           (3 * kW_BOSON_MASS * temp126 * temp28 * temp30 * temp8 *
            kW_BOSON_WIDTH) /
               4. +
           (3 * kW_BOSON_MASS * temp127 * temp28 * temp30 * temp8 *
            kW_BOSON_WIDTH) /
               4. +
           3 * kW_BOSON_MASS * temp133 * temp19 * temp5 * temp7 * temp8 *
               kW_BOSON_WIDTH -
           3 * kW_BOSON_MASS * temp134 * temp19 * temp5 * temp7 * temp8 *
               kW_BOSON_WIDTH +
           12 * kW_BOSON_MASS * temp123 * temp5 * kW_BOSON_WIDTH *
               std::atan(temp136 * temp67 * temp70) -
           12 * kW_BOSON_MASS * temp123 * temp5 * kW_BOSON_WIDTH *
               std::atan(temp138 * temp67 * temp70) +
           6 * kW_BOSON_MASS * temp5 * kW_BOSON_WIDTH *
               (temp5 * (-temp129 + temp1 * temp2 + temp132 * temp4 +
                         temp3 * (2 * temp132 + temp55) + temp35 * temp7 -
                         s * temp72 + temp1 * temp72 + 2 * s * temp1 * temp9 -
                         2 * temp2 * temp9 + s * temp73 * temp9 -
                         temp1 * temp73 * temp9 -
                         2 * kW_BOSON_MASS * kZ_BOSON_MASS * s *
                             kW_BOSON_WIDTH * kZ_BOSON_WIDTH -
                         2 * kZ_BOSON_MASS * temp71 * kW_BOSON_WIDTH *
                             kZ_BOSON_WIDTH) +
                temp17 *
                    (temp129 - temp1 * temp2 + s * temp72 - temp1 * temp72 -
                     2 * s * temp1 * temp9 + 2 * temp2 * temp9 -
                     s * temp73 * temp9 + temp1 * temp73 * temp9 +
                     temp35 * temp96 +
                     2 * kW_BOSON_MASS * kZ_BOSON_MASS * s * kW_BOSON_WIDTH *
                         kZ_BOSON_WIDTH +
                     2 * kZ_BOSON_MASS * temp71 * kW_BOSON_WIDTH *
                         kZ_BOSON_WIDTH +
                     temp4 * (s * temp7 + temp7 * temp9 -
                              kW_BOSON_MASS * kZ_BOSON_MASS * kW_BOSON_WIDTH *
                                  kZ_BOSON_WIDTH) +
                     temp3 * (temp55 + 2 * kW_BOSON_MASS *
                                           (-(kW_BOSON_MASS * s) +
                                            kW_BOSON_MASS * temp1 -
                                            kZ_BOSON_MASS * kW_BOSON_WIDTH *
                                                kZ_BOSON_WIDTH)))) *
               std::log(((temp128 + pow(temp136, 2)) *
                         (temp128 + (temp133 * temp8) / 4.)) /
                        ((temp128 + pow(temp138, 2)) *
                         (temp128 + (temp134 * temp8) / 4.))) +
           6 * kW_BOSON_MASS * temp139 * temp24 * temp65 * kW_BOSON_WIDTH *
               std::log(
                   (temp128 + temp141 + temp72 + temp34 * temp40 * temp9) /
                   (temp128 + temp140 + temp72 - temp23 * temp40 * temp9)) -
           6 * kW_BOSON_MASS * temp139 * temp24 * temp65 * kW_BOSON_WIDTH *
               std::log(
                   (temp140 + temp72 +
                    (s + temp12 - 2 * temp3 - temp16 * temp40 + temp73) *
                        temp9) /
                   (temp141 + temp72 +
                    temp40 * (temp16 + temp2 + temp31 + temp32 + s * temp73) *
                        temp9)))) /
         (24. * pow(kCOS_THETA_WEAK, 4) * pow(kSIN_THETA_WEAK, 4) * temp65);
}

/**
 * Compute the partially integrated squared matrix element for the decay of a
 * right-handed neutrino into an three active neutrino, all of the same
 * generation.
 * @param s Mandelstam variable s = (pu + pubar)^2.
 * @param params Contains mvr, theta.
 */
auto partialy_integrated_msqrd_vr_to_vl_vl_vl(double s, void *params)
    -> double {

  auto *pars = static_cast<PartiallyIntegratedMsqrdParams *>(params);
  const double mvr = pars->mvr;
  const double theta = pars->theta;

  double temp1 = pow(mvr, 2);
  double temp2 = -s;
  double temp3 = temp1 + temp2;
  double temp4 = pow(kZ_BOSON_MASS, 4);
  double temp5 = pow(kZ_BOSON_MASS, 2);
  double temp6 = pow(s, 2);
  double temp7 = -2 * s;
  double temp8 = pow(kZ_BOSON_WIDTH, 2);
  double temp9 = temp7 + temp8;
  double temp10 = temp5 * temp9;
  double temp11 = temp10 + temp4 + temp6;
  double temp12 = 1 / temp11;
  double temp13 = temp2 + temp5;
  double temp14 = -temp5;
  double temp15 = 2 * temp6;
  double temp16 = 2 * s;
  double temp17 = 1 / kZ_BOSON_WIDTH;
  double temp18 = -temp8;
  double temp19 = temp16 + temp18;
  double temp20 = temp19 * temp5;
  double temp21 = kZ_BOSON_MASS * temp17;
  double temp22 = std::atan(temp21);
  double temp23 = 1 / kZ_BOSON_MASS;
  double temp24 = temp16 + temp5;
  double temp25 = -(temp1 * temp24);
  double temp26 = temp15 + temp20 + temp25 + temp4;
  double temp27 = temp1 + temp14 + temp2;
  double temp28 = -2 * temp1 * temp5;
  double temp29 = 3 * temp4;
  double temp30 = -temp6;
  double temp31 = temp20 + temp28 + temp29 + temp30;
  double temp32 = -temp1;
  double temp33 = s + temp32 + temp5;
  double temp34 = -2 * temp6;
  double temp35 = temp17 * temp23 * temp33;
  double temp36 = std::atan(temp35);
  double temp37 = -(temp5 * temp8);
  double temp38 = pow(temp33, 2);
  double temp39 = temp5 * temp8;
  return (pow(kQE, 4) * pow(theta, 2) *
          (2 * temp17 * temp22 * temp23 * temp26 + 2 * temp3 -
           s * temp12 * pow(temp3, 2) - (temp12 * pow(temp3, 3)) / 3. -
           2 * temp12 * temp13 * temp38 +
           4 * temp12 * temp27 *
               (temp1 * temp13 + temp15 + temp39 - 2 * temp4) +
           temp17 * temp23 * temp36 *
               (temp10 + temp1 * temp24 + temp34 - temp4) +
           2 * temp12 * temp13 * temp4 -
           4 * temp12 * (temp1 * (s + temp14) + temp34 + temp37 + 2 * temp4) *
               temp5 +
           4 * kZ_BOSON_MASS * temp12 * temp22 * temp31 * kZ_BOSON_WIDTH -
           4 * kZ_BOSON_MASS * temp12 * temp31 * temp36 * kZ_BOSON_WIDTH +
           temp17 * temp23 * temp26 * std::atan(temp17 * temp23 * temp27) +
           (temp1 + (2 * s * temp3) / (temp1 + temp2 - 2 * temp5) -
            2 * (s + temp5) -
            2 * temp12 *
                (-pow(kZ_BOSON_MASS, 6) + pow(s, 3) +
                 temp1 * (temp30 + temp37 + temp4) - temp4 * (s - 3 * temp8) +
                 s * temp5 * (s + temp8))) *
               log((temp5 * (temp5 + temp8)) / (temp38 + temp39)))) /
         (4. * pow(kCOS_THETA_WEAK, 4) * pow(kSIN_THETA_WEAK, 4));
}

} // namespace storm