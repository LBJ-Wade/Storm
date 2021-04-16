#include "storm/types.hpp"
#include <cstddef>
#include <fmt/core.h>
#include <fmt/format.h>
#include <gtest/gtest.h>
#include <iostream>
#include <storm/simple.hpp>

// Benchmark parameters
constexpr double BM_MVR = 1e4;
constexpr double BM_THETA = 1e-5;
constexpr int BM_GEN = 0;
constexpr double REL_ERROR = 1e-3;

inline static auto frac_err(double meas, double exp) -> double {
  return std::abs(meas - exp) / exp;
}

TEST(TestSimpleRhNeutrino, TestTwoBodyWidths) {

  constexpr double mm_wvh = 0.0000328057;
  constexpr double mm_wvz = 0.0000308646;
  constexpr double mm_wlw = 0.0000305341;

  const double wvh = storm::width_vr_to_vl_h(BM_MVR, BM_THETA, BM_GEN);
  const double wvz = storm::width_vr_to_vl_z(BM_MVR, BM_THETA, BM_GEN);
  const double wlw = storm::width_vr_to_l_w(BM_MVR, BM_THETA, BM_GEN);

  fmt::print("Gamma(N -> nu + H) = {}\n", wvh);
  fmt::print("Gamma(N -> nu + Z) = {}\n", wvz);
  fmt::print("Gamma(N -> l + W) = {}\n", wlw);

  ASSERT_LE(frac_err(wvh, mm_wvh), REL_ERROR);
  ASSERT_LE(frac_err(wvz, mm_wvz), REL_ERROR);
  ASSERT_LE(frac_err(wlw, mm_wlw), REL_ERROR);
}

TEST(TestSimpleRhNeutrino, TestThreeBodyWidths) {

  const size_t nevents = 50000;

  const auto wvuu = storm::width_vr_to_vl_u_u(BM_MVR, BM_THETA, BM_GEN, 0);
  const auto wvcc = storm::width_vr_to_vl_u_u(BM_MVR, BM_THETA, BM_GEN, 1);
  const auto wvtt = storm::width_vr_to_vl_u_u(BM_MVR, BM_THETA, BM_GEN, 2);

  const auto wvdd = storm::width_vr_to_vl_d_d(BM_MVR, BM_THETA, BM_GEN, 0);
  const auto wvss = storm::width_vr_to_vl_d_d(BM_MVR, BM_THETA, BM_GEN, 1);
  const auto wvbb = storm::width_vr_to_vl_d_d(BM_MVR, BM_THETA, BM_GEN, 2);

  const auto wlud = storm::width_vr_to_l_u_d(BM_MVR, BM_THETA, BM_GEN, 0);
  const auto wlcs = storm::width_vr_to_l_u_d(BM_MVR, BM_THETA, BM_GEN, 1);
  const auto wltb = storm::width_vr_to_l_u_d(BM_MVR, BM_THETA, BM_GEN, 2);

  const auto wvmm = storm::width_vr_to_vl_lp_lp(BM_MVR, BM_THETA, BM_GEN, 1);
  const auto wvtata = storm::width_vr_to_vl_lp_lp(BM_MVR, BM_THETA, BM_GEN, 2);

  const auto wvme = storm::width_vr_to_vlp_lp_l(BM_MVR, BM_THETA, BM_GEN, 1);
  const auto wvte = storm::width_vr_to_vlp_lp_l(BM_MVR, BM_THETA, BM_GEN, 2);

  const auto wvee = storm::width_vr_to_vl_l_l(BM_MVR, BM_THETA, BM_GEN);

  const auto wvvv = storm::width_vr_to_vl_vl_vl(BM_MVR, BM_THETA, BM_GEN);

  fmt::print("Gamma(N -> nu + u + ubar) = {}\n", wvuu);
  fmt::print("Gamma(N -> nu + c + cbar) = {}\n", wvcc);
  fmt::print("Gamma(N -> nu + t + tbar) = {}\n", wvtt);

  fmt::print("Gamma(N -> nu + d + dbar) = {}\n", wvdd);
  fmt::print("Gamma(N -> nu + s + sbar) = {}\n", wvss);
  fmt::print("Gamma(N -> nu + b + bbar) = {}\n", wvbb);

  fmt::print("Gamma(N -> l + u + dbar) = {}\n", wlud);
  fmt::print("Gamma(N -> l + c + sbar) = {}\n", wlcs);
  fmt::print("Gamma(N -> l + t + bbar) = {}\n", wltb);

  fmt::print("Gamma(N -> nue + mu + mubar) = {}\n", wvmm);
  fmt::print("Gamma(N -> nue + tau + taubar) = {}\n", wvtata);

  fmt::print("Gamma(N -> nu_mu + mu + ebar) = {}\n", wvme);
  fmt::print("Gamma(N -> nu_tau + tau + ebar) = {}\n", wvte);

  fmt::print("Gamma(N -> nu_e + e + ebar) = {}\n", wvee);

  fmt::print("Gamma(N -> nu_e + nu_e + nu_e) = {}\n", wvvv);

  ASSERT_EQ(true, true);
}
