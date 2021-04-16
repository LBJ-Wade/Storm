"""
Module containing the simpler definitions of the spectra functions for the
`SimpleRhNeutrino` model.
"""
import numpy as np
from storm.models.simple._simple_spectra import (
    spectrum_vr_to_vl_h as _dndx_vl_h,
    spectrum_vr_to_vl_z as _dndx_vl_z,
    spectrum_vr_to_l_w as _dndx_l_w,
    spectrum_vr_to_vl_u_u as _dndx_vl_u_u,
    spectrum_vr_to_vl_d_d as _dndx_vl_d_d,
    spectrum_vr_to_l_u_d as _dndx_l_u_d,
    spectrum_vr_to_vl_l_l as _dndx_vl_l_l,
    spectrum_vr_to_vl_lp_lp as _dndx_vl_lp_lp,
    spectrum_vr_to_vlp_lp_l as _dndx_vlp_lp_l,
)


def dndx_vl_h(**kwargs):
    spec = _dndx_vl_h(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
    )
    return np.array(spec[0]), np.array(spec[1])


def dndx_vl_z(**kwargs):
    spec = _dndx_vl_z(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
    )
    return np.array(spec[0]), np.array(spec[1])


def dndx_l_w(**kwargs):
    spec = _dndx_l_w(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
        anti=kwargs["anti"],
    )
    return np.array(spec[0]), np.array(spec[1])


def dndx_vl_u_u(**kwargs):
    spec = _dndx_vl_u_u(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        genq=kwargs["genq"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
    )
    return np.array(spec[0]), np.array(spec[1])


def dndx_vl_d_d(**kwargs):
    spec = _dndx_vl_d_d(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        genq=kwargs["genq"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
    )
    return np.array(spec[0]), np.array(spec[1])


def dndx_l_u_d(**kwargs):
    spec = _dndx_l_u_d(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        genq=kwargs["genq"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
        anti=kwargs["anti"],
    )
    return np.array(spec[0]), np.array(spec[1])


def dndx_vl_l_l(**kwargs):
    spec = _dndx_vl_l_l(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
    )
    return np.array(spec[0]), np.array(spec[1])


def dndx_vl_lp_lp(**kwargs):
    spec = _dndx_vl_lp_lp(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        genlp=kwargs["genlp"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
    )
    return np.array(spec[0]), np.array(spec[1])


def dndx_vlp_lp_l(**kwargs):
    spec = _dndx_vlp_lp_l(
        mvr=kwargs["mvr"],
        theta=kwargs["theta"],
        genl=kwargs["genl"],
        genlp=kwargs["genlp"],
        product=kwargs["product"],
        xbounds=kwargs["xbounds"],
        nbins=kwargs["nbins"],
        nevents=kwargs["nevents"],
        anti=kwargs["anti"],
    )
    return np.array(spec[0]), np.array(spec[1])
