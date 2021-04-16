"""
Module containing the definitions of all the decay widths for the 
`SimpleRhNeutrino` model.
"""

from storm.models.simple._simple_widths import (
    width_vr_to_l_w as _width_vr_to_l_w,
    width_vr_to_vl_h as _width_vr_to_vl_h,
    width_vr_to_vl_z as _width_vr_to_vl_z,
    width_vr_to_vl_u_u as _width_vr_to_vl_u_u,
    width_vr_to_vl_d_d as _width_vr_to_vl_d_d,
    width_vr_to_l_u_d as _width_vr_to_l_u_d,
    width_vr_to_vl_lp_lp as _width_vr_to_vl_lp_lp,
    width_vr_to_vlp_lp_l as _width_vr_to_vlp_lp_l,
    width_vr_to_vl_l_l as _width_vr_to_vl_l_l,
    width_vr_to_vl_vl_vl as _width_vr_to_vl_vl_vl,
)


def width_vl_h(mvr, theta, genl):
    return _width_vr_to_vl_h(mvr=mvr, theta=theta, genl=genl)


def width_vl_z(mvr, theta, genl):
    return _width_vr_to_vl_z(mvr=mvr, theta=theta, genl=genl)


def width_l_w(mvr, theta, genl):
    return _width_vr_to_l_w(mvr=mvr, theta=theta, genl=genl)


def width_vl_u_u(mvr, theta, genl, genq):
    return _width_vr_to_vl_u_u(mvr=mvr, theta=theta, genl=genl, genq=genq)


def width_vl_d_d(mvr, theta, genl, genq):
    return _width_vr_to_vl_d_d(mvr=mvr, theta=theta, genl=genl, genq=genq)


def width_l_u_d(mvr, theta, genl, genq):
    return _width_vr_to_l_u_d(mvr=mvr, theta=theta, genl=genl, genq=genq)


def width_vl_lp_lp(mvr, theta, genl, genlp):
    return _width_vr_to_vl_lp_lp(mvr=mvr, theta=theta, genl=genl, genlp=genlp)


def width_vlp_lp_l(mvr, theta, genl, genlp):
    return _width_vr_to_vlp_lp_l(mvr=mvr, theta=theta, genl=genl, genlp=genlp)


def width_vl_l_l(mvr, theta, genl):
    return _width_vr_to_vl_l_l(mvr=mvr, theta=theta, genl=genl)


def width_vl_vl_vl(mvr, theta, genl):
    return _width_vr_to_vl_vl_vl(mvr=mvr, theta=theta, genl=genl)
