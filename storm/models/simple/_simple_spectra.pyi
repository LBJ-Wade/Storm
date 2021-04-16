from typing import Tuple, List

def spectrum_vr_to_vl_h(
    mvr: float,
    theta: float,
    genl: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
) -> Tuple[List[float], List[float]]: ...
def spectrum_vr_to_vl_z(
    mvr: float,
    theta: float,
    genl: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
) -> Tuple[List[float], List[float]]: ...
def spectrum_vr_to_l_w(
    mvr: float,
    theta: float,
    genl: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
    anti: bool,
) -> Tuple[List[float], List[float]]: ...
def spectrum_vr_to_vl_u_u(
    mvr: float,
    theta: float,
    genl: int,
    genq: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
) -> Tuple[List[float], List[float]]: ...
def spectrum_vr_to_vl_d_d(
    mvr: float,
    theta: float,
    genl: int,
    genq: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
) -> Tuple[List[float], List[float]]: ...
def spectrum_vr_to_l_u_d(
    mvr: float,
    theta: float,
    genl: int,
    genq: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
    anti: bool,
) -> Tuple[List[float], List[float]]: ...
def spectrum_vr_to_vl_l_l(
    mvr: float,
    theta: float,
    genl: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
) -> Tuple[List[float], List[float]]: ...
def spectrum_vr_to_vl_lp_lp(
    mvr: float,
    theta: float,
    genl: int,
    genlp: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
) -> Tuple[List[float], List[float]]: ...
def spectrum_vr_to_vlp_lp_l(
    mvr: float,
    theta: float,
    genl: int,
    genlp: int,
    product: int,
    xbounds: Tuple[float, float],
    nbins: int,
    nevents: int,
    anti: bool,
) -> Tuple[List[float], List[float]]: ...