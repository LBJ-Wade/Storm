#include "storm/simple.hpp"
#include <fmt/format.h>
#include <gtest/gtest.h>
#include <storm/constants.hpp>
#include <storm/pythia_spectra.hpp>

constexpr double BM_MVR = 1e4;
constexpr double BM_THETA = 1e-4;
constexpr int BM_GEN = 0;

TEST(TestSpectra, TestTwoBody) {
  const int product = storm::kPHOTON_ID;
  const auto xbounds = std::make_pair(1e-4, 1.0);
  const unsigned int nbins = 150;
  const size_t nevents = 50000;

  auto spectrum_vl_h = storm::spectrum_vr_to_vl_h(
      BM_MVR, BM_THETA, BM_GEN, product, xbounds, nbins, nevents);
  auto &xs = spectrum_vl_h.first;
  auto &spec = spectrum_vl_h.second;

  for (size_t i = 0; i < xs.size(); i++) {
    const double x = xs[i];
    const double dndx = spec[i];

    auto res = fmt::format("{}, {}", x, dndx);
    std::cout << "{" + res + "},\n";
  }
}