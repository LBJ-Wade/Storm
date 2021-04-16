"""
Module containing various physical constants.
"""

import numpy as np

G_FERMI = 1.1663787e-5
HIGGS_VEV = 246.21965
ALPHA_EM = 1.0 / 137.0
# at p^2 = 0
QE = 0.302862
SIN_THETA_WEAK = 0.480853
SIN_THETA_WEAK_SQRD = 0.23122
COS_THETA_WEAK = 0.876801
M_PLANK = 1.220910e19

# Masses
ELECTRON_MASS = 0.5109989461e-3
MUON_MASS = 105.6583745e-3
TAU_MASS = 1776.86e-3

LEPTON_MASSES = np.array([ELECTRON_MASS, MUON_MASS, TAU_MASS], dtype=np.float_)

UP_QUARK_MASS = 2.16e-3
CHARM_QUARK_MASS = 1.27
TOP_QUARK_MASS = 172.9

UP_QUARK_MASSES = np.array(
    [UP_QUARK_MASS, CHARM_QUARK_MASS, TOP_QUARK_MASS], dtype=np.float_
)

DOWN_QUARK_MASS = 4.67e-3
STRANGE_QUARK_MASS = 93e-3
BOTTOM_QUARK_MASS = 4.18

DOWN_QUARK_MASSES = np.array(
    [DOWN_QUARK_MASS, STRANGE_QUARK_MASS, BOTTOM_QUARK_MASS], dtype=np.float_
)

W_BOSON_MASS = 80.38500
Z_BOSON_MASS = 91.18760
HIGGS_MASS = 125.00
NEUTRAL_PION_MASS = 134.9766e-3
CHARGED_PION_MASS = 139.57018e-3
NEUTRAL_KAON_MASS = 497.61e-3
LONG_KAON_MASS = 497.614e-3
CHARGED_KAON_MASS = 493.68e-3

# Widths
W_BOSON_WIDTH = 2.08500
Z_BOSON_WIDTH = 2.49520
HIGGS_WIDTH = 0.00374

# CKM
VCKM = np.array(
    [
        [
            0.9742855469766043 + 0j,
            0.22528978739658034 + 0j,
            0.003467579534847188 - 0.000400673772678475j,
        ],
        [
            -0.22523919033512343 + 0.000016074829716163196j,
            0.9734329404782395 + 3.7170775861643788e-6j,
            0.041177873389110276 + 0.0j,
        ],
        [
            0.005901499687662993 + 0.0003900419379730849j,
            -0.04090004814585696 + 0.0000901916953960691j,
            0.9991457341628637,
        ],
    ],
    dtype=np.complex_,
)