#ifndef STORM_CONSTANTS_HPP
#define STORM_CONSTANTS_HPP

#include <array>
#include <complex>

namespace storm {
// ==========================
// ---- Physical Constants --
// ==========================

// Various physical constants
static constexpr double kG_FERMI = 1.1663787e-5;
static constexpr double kHIGGS_VEV = 246.21965;
static constexpr double kALPHA_EM = 1.0 / 137.0; // at p^2 = 0
static constexpr double kQE = 0.302862;
static constexpr double kSIN_THETA_WEAK = 0.480853;
static constexpr double kSIN_THETA_WEAK_SQRD = 0.23122;
static constexpr double kCOS_THETA_WEAK = 0.876801;
static constexpr double kM_PLANK = 1.220910e19;
static constexpr double kRHO_CRIT = 1.05375e-5;
static constexpr double kS_TODAY = 2891.2;
static constexpr double kT_CMB = 2.56215e-10;
static constexpr double kT_BBN = 0.0001; // 0.1 MeV in GeV
static constexpr double kOMEGA_H2_CDM = 0.1198;

// Masses
static constexpr double kELECTRON_MASS = 0.5109989461e-3;
static constexpr double kMUON_MASS = 105.6583745e-3;
static constexpr double kTAU_MASS = 1776.86e-3;
static constexpr double kUP_QUARK_MASS = 2.16e-3;
static constexpr double kDOWN_QUARK_MASS = 4.67e-3;
static constexpr double kSTRANGE_QUARK_MASS = 93e-3;
static constexpr double kCHARM_QUARK_MASS = 1.27;
static constexpr double kBOTTOM_QUARK_MASS = 4.18;
static constexpr double kTOP_QUARK_MASS = 172.9;
static constexpr double kW_BOSON_MASS = 80.38500;
static constexpr double kZ_BOSON_MASS = 91.18760;
static constexpr double kHIGGS_MASS = 125.00;
static constexpr double kNEUTRAL_PION_MASS = 134.9766e-3;
static constexpr double kCHARGED_PION_MASS = 139.57018e-3;
static constexpr double kNEUTRAL_KAON_MASS = 497.61e-3;
static constexpr double kLONG_KAON_MASS = 497.614e-3;
static constexpr double kCHARGED_KAON_MASS = 493.68e-3;

static constexpr std::array<double, 3> kLEPTON_MASSES = {kELECTRON_MASS,
                                                         kMUON_MASS, kTAU_MASS};
static constexpr std::array<double, 3> kUP_QUARK_MASSES = {
    kUP_QUARK_MASS, kCHARM_QUARK_MASS, kTOP_QUARK_MASS};
static constexpr std::array<double, 3> kDOWN_QUARK_MASSES = {
    kDOWN_QUARK_MASS, kSTRANGE_QUARK_MASS, kBOTTOM_QUARK_MASS};

// Boson widths
static constexpr double kW_BOSON_WIDTH = 2.08500;
static constexpr double kZ_BOSON_WIDTH = 2.49520;
static constexpr double kHIGGS_WIDTH = 0.00374;

// PDG Codes
static constexpr int kELECTRON_ID = 11;
static constexpr int kELECTRON_NEUTRINO_ID = 12;
static constexpr int kMUON_ID = 13;
static constexpr int kMUON_NEUTRINO_ID = 14;
static constexpr int kTAU_ID = 15;
static constexpr int kTAU_NEUTRINO_ID = 16;
static constexpr int kHIGGS_ID = 25;
static constexpr int kZ_BOSON_ID = 23;
static constexpr int kW_BOSON_ID = 24;
static constexpr int kPHOTON_ID = 22;
static constexpr int kUP_QUARK_ID = 2;
static constexpr int kCHARM_QUARK_ID = 4;
static constexpr int kTOP_QUARK_ID = 6;
static constexpr int kDOWN_QUARK_ID = 1;
static constexpr int kSTRANGE_QUARK_ID = 3;
static constexpr int kBOTTOM_QUARK_ID = 5;

// CKM Matrix
static constexpr std::array<std::array<std::complex<double>, 3>, 3> kCKM = {
    std::array<std::complex<double>, 3>{
        std::complex<double>{0.9742855469766043, 0.0},
        std::complex<double>{0.22528978739658034, 0.0},
        std::complex<double>{0.003467579534847188, -0.000400673772678475}},
    std::array<std::complex<double>, 3>{
        std::complex<double>{-0.22523919033512343, 0.000016074829716163196},
        std::complex<double>{0.9734329404782395, 3.7170775861643788e-6},
        std::complex<double>{0.041177873389110276, 0.0}},
    std::array<std::complex<double>, 3>{
        std::complex<double>{0.005901499687662993, 0.0003900419379730849},
        std::complex<double>{-0.04090004814585696, 0.0000901916953960691},
        std::complex<double>{0.9991457341628637}},
};

} // namespace storm

#endif // STORM_CONSTANTS_HPP
