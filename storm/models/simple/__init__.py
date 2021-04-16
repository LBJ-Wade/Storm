"""
Module for describing a RH-Neutrino which couples only with a single active
neutrino through a Yukawa interaction with the SM Higgs.
"""
from abc import ABC, abstractmethod
from functools import partial
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np
from hazma.rh_neutrino import RHNeutrino as _HazmaRhNeutrino  # type: ignore

from storm.constants import LEPTON_MASSES
from storm.models.simple._spectra import dndx_l_u_d as _dndx_l_u_d
from storm.models.simple._spectra import dndx_l_w as _dndx_l_w
from storm.models.simple._spectra import dndx_vl_d_d as _dndx_vl_d_d
from storm.models.simple._spectra import dndx_vl_h as _dndx_vl_h
from storm.models.simple._spectra import dndx_vl_l_l as _dndx_vl_l_l
from storm.models.simple._spectra import dndx_vl_lp_lp as _dndx_vl_lp_lp
from storm.models.simple._spectra import dndx_vl_u_u as _dndx_vl_u_u
from storm.models.simple._spectra import dndx_vl_z as _dndx_vl_z
from storm.models.simple._spectra import dndx_vlp_lp_l as _dndx_vlp_lp_l
from storm.models.simple._widths import width_l_u_d as _width_l_u_d
from storm.models.simple._widths import width_l_w as _width_l_w
from storm.models.simple._widths import width_vl_d_d as _width_vl_d_d
from storm.models.simple._widths import width_vl_h as _width_vl_h
from storm.models.simple._widths import width_vl_l_l as _width_vl_l_l
from storm.models.simple._widths import width_vl_lp_lp as _width_vl_lp_lp
from storm.models.simple._widths import width_vl_u_u as _width_vl_u_u
from storm.models.simple._widths import width_vl_vl_vl as _width_vl_vl_vl
from storm.models.simple._widths import width_vl_z as _width_vl_z
from storm.models.simple._widths import width_vlp_lp_l as _width_vlp_lp_l

# Type of the final states: ('s1', 's2', ...)
StateType = Tuple[str, ...]

_PRODUCT_NAME_TO_ID = {
    "photon": 22,
    "electron": 11,
    "positron": -11,
    "electron-neutrino": 12,
    "muon-neutrino": 14,
    "tau-neutrino": 16,
    "neutron": 2112,
    "proton": 2212,
    "anti-proton": -2212,
}

AVAILIBLE_PRODUCTS = _PRODUCT_NAME_TO_ID.keys()

# Type of the underlying c++ functions.
_SPECTRUM_CALLABLE = Callable[
    [float, float, int, int, Tuple[float, float], int, int],
    Tuple[List[float], List[float]],
]


class SimpleRhNeutrinoPythia(SimpleRhNeutrinoBase):
    """"""

    def __init__(self, mvr: float, theta: float, lep: str):
        super().__init__(mvr, theta, lep)

        # Dictionary of functions to compute the partial widths of a given
        # final state.
        self._width_dispatch: Dict[Tuple[str, ...],
                                   Callable[..., float]] = dict()
        # Dictionary of functions to compute the spectra of from the decay into
        # a given final state.
        self._dndx_dispatch: Dict[
            Tuple[str, ...],
            Callable[..., Tuple[np.ndarray, np.ndarray]]] = dict()
        # List of tuples specifying all decay modes
        self._decay_final_states: List[Tuple[str, ...]] = list()
        # Dictionary specifying to conjugate of a given final state.
        self._conj_map: Dict[Tuple[str, ...], Tuple[str, ...]] = dict()

        self.__create_states_and_dispatch_tables()

    def __create_states_and_dispatch_tables(self):
        """
        Create the dispatch tables for the partial width and spectrum
        functions.
        """
        lep = self._lep
        genl = self._genl
        dt = {
            (f"v{lep}", "h"): partial(_width_vl_h, genl=genl),
            (f"v{lep}", "z"): partial(_width_vl_z, genl=genl),
            (f"{lep}", "w"): partial(_width_l_w, genl=genl),
            (f"{lep}bar", "wbar"): partial(_width_l_w, genl=genl),
            (f"v{lep}", f"v{lep}", f"v{lep}"):
                partial(_width_vl_vl_vl, genl=genl),
            (f"v{lep}", f"{lep}", f"{lep}bar"):
                partial(_width_vl_l_l, genl=genl),
        }
        sdt = {
            (f"v{lep}", "h"): partial(_dndx_vl_h, genl=genl),
            (f"v{lep}", "z"): partial(_dndx_vl_z, genl=genl),
            (f"{lep}", "w"): partial(_dndx_l_w, genl=genl, anti=False),
            (f"{lep}bar", "wbar"): partial(_dndx_l_w, genl=genl, anti=True),
            (f"v{lep}", f"{lep}", f"{lep}bar"):
                partial(_dndx_vl_l_l, genl=genl),
        }

        self._conj_map = {(f"{lep}bar", "wbar"): (f"{lep}", "w")}

        # Add states of the form vl + q + qbar and l + u + dbar
        for i, (u, d) in enumerate([("u", "d"), ("c", "s"), ("t", "b")]):
            state = (f"v{lep}", u, f"{u}bar")
            dt[state] = partial(_width_vl_u_u, genl=genl, genq=i)
            sdt[state] = partial(_dndx_vl_u_u, genl=genl, genq=i)

            state = (f"v{lep}", d, f"{d}bar")
            dt[state] = partial(_width_vl_d_d, genl=genl, genq=i)
            sdt[state] = partial(_dndx_vl_d_d, genl=genl, genq=i)

            state = (f"{lep}", u, f"{d}bar")
            dt[state] = partial(_width_l_u_d, genl=genl, genq=i)
            sdt[state] = partial(_dndx_l_u_d, genl=genl, genq=i, anti=False)

            state = (f"{lep}bar", f"{u}bar", f"{d}")
            dt[state] = partial(_width_l_u_d, genl=genl, genq=i)
            sdt[state] = partial(_dndx_l_u_d, genl=genl, genq=i, anti=True)

            self._conj_map[state] = (f"{lep}", u, f"{d}bar")

        for i, ell in enumerate(["e", "mu", "tau"]):
            if ell != lep:
                state = (f"v{lep}", f"{ell}", f"{ell}bar")
                dt[state] = partial(_width_vl_lp_lp, genl=genl, genlp=i)
                sdt[state] = partial(_dndx_vl_lp_lp, genl=genl, genlp=i)

                state = (f"v{ell}", f"{ell}", f"{lep}bar")
                dt[state] = partial(_width_vlp_lp_l, genl=genl, genlp=i)
                sdt[state] = partial(
                    _dndx_vlp_lp_l, genl=genl, genlp=i, anti=False)

                state = (f"v{ell}", f"{lep}", f"{ell}bar")
                dt[state] = partial(_width_vlp_lp_l, genl=genl, genlp=i)
                sdt[state] = partial(
                    _dndx_vlp_lp_l, genl=genl, genlp=i, anti=True)

                self._conj_map[state] = (f"v{ell}", f"{ell}", f"{lep}bar")

        self._width_dispatch = dt
        self._dndx_dispatch = sdt
        self._decay_final_states = list(dt.keys())

    @property
    def lep(self) -> str:
        return self._lep

    @lep.setter
    def lep(self, lep: str) -> None:
        self._lep = lep
        self.__create_states_and_dispatch_tables()

    @property
    def decay_final_states(self) -> List[StateType]:
        return self._decay_final_states

    @decay_final_states.setter
    def decay_final_states(self, val: List[StateType]) -> None:
        raise AttributeError('Cannot set "decay final states."')

    def partial_width(self, state: StateType, **kwargs) -> float:
        """
        Compute the partial width for a right-handed neutrino to decay into a
        particular final state.

        Parameters
        ----------
        state: Tuple[str, ...]
            A tuple of strings representing the final state. For example,
            state=('ve', 'h') will return the partial width for vr -> ve + h.
            See `decay_final_states` for a list of all final states.

        Returns
        -------
        pw: float
            The partial decay width.
        """
        if state in self._width_dispatch.keys():
            return self._width_dispatch[state](self._mvr, self._theta)
        raise ValueError(f"Invalid state: {state}")

    def partial_widths(self, **kwargs) -> Dict[StateType, float]:
        """
        Compute the partial width for a right-handed neutrino to decay into all
        possible final states.

        Parameters
        ----------
        remove_conjugates: Optional[bool]
            If true, the conjugate states are removed and their partial widths
            are added into the unconjugated state.

        Returns
        -------
        pws: Dict[str, float]
            Dictionary containing all partial decay widths from the decay of a
            right-handed neutrino.
        """
        pws = {
            key: func(self._mvr, self._theta)
            for key, func in self._width_dispatch.items()
        }

        if "remove_conjugates" in kwargs:
            for key, val in self._conj_map.items():
                pws[val] += pws[key]
                del pws[key]

        pws[("total",)] = sum(pws.values())
        return pws

    def branching_fractions(self, **kwargs) -> Dict[StateType, float]:
        """
        Compute the branching fractions for a right-handed neutrino to decay
        into a all availible final states.

        Parameters
        ----------
        remove_conjugates: Optional[bool]
            If true, the conjugate states are removed and their partial widths
            are added into the unconjugated state.

        Returns
        -------
        bf: Dict[str, float]
            Dictionary containing all branching fractions from the decay of a
            right-handed neutrino.
        """
        remove_conjugates = kwargs.get("remove_conjugates")
        if remove_conjugates is None:
            remove_conjugates = True

        pws = self.partial_widths(remove_conjugates=remove_conjugates)
        return {
            key: val / pws[("total",)]
            for key, val in pws.items() if key != ("total",)
        }

    def dndx_single_state(
        self,
        x: np.ndarray,
        product: int,
        state: StateType,
        **kwargs
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the spectrum of a specified product from the decay of a
        right-handed neutrino into all availible final states of a specified
        final state.

        Parameters
        ----------
        product: int,
            PDG code of the product to compute spectrum for. For example, to
            compute the photon spectrum, use `22`.
        xbounds: Tuple[float, float]
            Bounds on `x = 2*E/mvr`.
        nevents: Optional[int]
            Number of Pythia events to use in generating the spectrum. Default
            is 10_000.
        state: Optional[Tuple[str, ...]]
            A tuple of strings representing the final state. For example,
            state=('ve', 'h') will return the partial width for vr -> ve + h.
            See `decay_final_states` for a list of all final states.

        Returns
        -------
        dndx:
            If a state was specified, the return is the x and spectrum values.
            Otherwise, the x values are returns along with a dictionary of the
            spectra for all posible final states.
        """

        nevents = kwargs["nevents"] if "nevents" in kwargs else 10_000

        # Set the arguments to pass to c++ functions
        _kwargs = {
            "mvr": self._mvr,
            "theta": self._theta,
            "product": product,
            "xbounds": (np.min(x), np.max(x)),
            "nbins": len(x),
            "nevents": nevents,
        }

        if state in self._dndx_dispatch.keys():
            return self._dndx_dispatch[state](**_kwargs)
        raise ValueError(f"Invalid state: {state}")

    def dndx(
            self,
            x: np.ndarray,
            product: int,
            **kwargs
    ) -> Tuple[np.ndarray, Dict[StateType, np.ndarray]]:
        """
        Compute the spectrum of a specified product from the decay of a
        right-handed neutrino into all availible final states of a specified
        final state.

        Parameters
        ----------
        product: int,
            PDG code of the product to compute spectrum for. For example, to
            compute the photon spectrum, use `22`.
        xbounds: Tuple[float, float]
            Bounds on `x = 2*E/mvr`.
        nevents: Optional[int]
            Number of Pythia events to use in generating the spectrum. Default
            is 10_000.
        state: Optional[Tuple[str, ...]]
            A tuple of strings representing the final state. For example,
            state=('ve', 'h') will return the partial width for vr -> ve + h.
            See `decay_final_states` for a list of all final states.

        Returns
        -------
        dndx:
            If a state was specified, the return is the x and spectrum values.
            Otherwise, the x values are returns along with a dictionary of the
            spectra for all posible final states.
        """

        nevents = kwargs["nevents"] if "nevents" in kwargs else 10_000

        # Set the arguments to pass to c++ functions
        _kwargs = {
            "mvr": self._mvr,
            "theta": self._theta,
            "product": product,
            "xbounds": (np.min(x), np.max(x)),
            "nbins": len(x),
            "nevents": nevents,
        }

        dndx = {key: np.zeros_like(x) for key in self._dndx_dispatch.keys()}
        # Use the first state to get xs
        first_state = list(self._dndx_dispatch.keys())[0]
        xs, dndx[first_state] = self._dndx_dispatch[first_state](**_kwargs)

        # Compute spectra for all other states
        dndx = {
            key: func(**_kwargs)[1]
            for key, func in self._dndx_dispatch.items()
            if key != first_state
        }

        # Apply branching fractions and compute total spectrum
        total = np.zeros_like(xs)
        bfs = self.branching_fractions(remove_conjugates=False)
        for key in dndx.keys():
            dndx[key] *= bfs[key]
            total += dndx[key]
        dndx[("total",)] = total

        return xs, dndx


class SimpleRhNeutrinoHazma(SimpleRhNeutrinoBase):
    """
    Model of a right-handed neutrino with a mass less than 1 GeV that mixes
    with a single active neutrino.
    """

    def __init__(self, mvr, theta, lep):
        super().__init__(mvr, theta, lep)

        self._hazma = _HazmaRhNeutrino(mvr, theta, lep, include_3body=True)

        # Dictionary of functions to compute the partial widths of a given
        # final state.
        self._width_dispatch: Dict[Tuple[str, ...],
                                   Callable[..., float]] = dict()
        # Dictionary of functions to compute the spectra of from the decay into
        # a given final state.
        self._dnde_dispatch: Dict[
            Tuple[str, ...],
            Callable[..., Tuple[np.ndarray, np.ndarray]]] = dict()
        # List of tuples specifying all decay modes
        self._decay_final_states: List[Tuple[str, ...]] = list()
        # Dictionary specifying to conjugate of a given final state.
        self._conj_map: Dict[Tuple[str, ...], Tuple[str, ...]] = dict()

    def __create_states_and_dispatch_tables(self):
        lep = self._lep

        self._width_dispatch = {
            (f"{lep}", "pi"): self._hazma.width_pi_l,
            (f"{lep}", "k"): self._hazma.width_k_l(),
            (f"v{lep}", "pi0"): self._hazma.width_pi0_nu,
            (f"v{lep}", "g"): self._hazma.width_nu_gamma,
            (f"v{lep}", "pi", "pi"): self._hazma.width_nu_pi_pi,
            (f"{lep}", "pi", "pi0"): self._hazma.width_l_pi_pi0(),
            (f"v{lep}", f"{lep}", f"{lep}"): self._hazma.width_nu_l_l,
            (f"v{lep}", f"v{lep}", f"v{lep}"): self._hazma.width_nu_nu_nu,
            (f"v{lep}", "g", "g"): self._hazma.width_nu_g_g,
        }

        self._dnde_dispatch = {
            (f"{lep}", "pi"): self._hazma.dnde_pi_l,
            (f"{lep}", "k"): self._hazma.dnde_k_l,
            (f"v{lep}", "pi0"): self._hazma.dnde_nu_pi0,
            # (f"v{lep}", "g"): self._hazma.dnde_,
            (f"v{lep}", "pi", "pi"): self._hazma.dnde_nu_pi_pi,
            (f"{lep}", "pi", "pi0"): self._hazma.dnde_l_pi_pi0,
            (f"v{lep}", f"{lep}", f"{lep}"): self._hazma.dnde_nu_l_l,
            (f"v{lep}", "g", "g"): self._hazma.dnde_nu_g_g,
        }

        self._decay_final_states = list(self._width_dispatch.keys())

        self._conj_map = {
            (f"{lep}", "pi"): (f"{lep}bar", "pibar"),
            (f"{lep}", "k"): (f"{lep}bar", "kbar"),
            (f"{lep}", "pi", "pi0"): (f"{lep}bar", "pibar", "pi0"),
        }

    def partial_width(self, state: StateType, **kwargs) -> float:
        if state in self._decay_final_states:
            return self._width_dispatch[state]()
        else:
            raise ValueError(f"Invalid state: {state}")

    def partial_widths(self, **kwargs) -> Dict[StateType, float]:
        widths = dict()
        for key, func in self._width_dispatch.items():
            if key not in self._conj_map:
                widths[key] = func()

        remove_conjugate = kwargs.get('remove_conjugate')
        if remove_conjugate is None:
            if not remove_conjugate:
                for state, conj in self._conj_map.items():
                    widths[state] /= 2.0
                    widths[conj] = widths[state]

        return {key: func() for key, func in self._width_dispatch.items()}

    def dndx(self, x: np.ndarray, product: int, **kwargs) -> float:
        """
        Compute the spectrum of a specified product from the decay of a
        right-handed neutrino into all availible final states of a specified
        final state.

        Parameters
        ----------
        product: int,
            PDG code of the product to compute spectrum for. For example, to
            compute the photon spectrum, use `22`.
        xbounds: Tuple[float, float]
            Bounds on `x = 2*E/mvr`.
        state: Optional[Tuple[str, ...]]
            A tuple of strings representing the final state. For example,
            state=('ve', 'h') will return the partial width for vr -> ve + h.
            See `decay_final_states` for a list of all final states.

        Returns
        -------
        dndx:
            If a state was specified, the return is the x and spectrum values.
            Otherwise, the x values are returns along with a dictionary of the
            spectra for all posible final states.
        """
        if not product == "photon" or not product == 22:
            raise NotImplementedError(
                "Only the photon spectrum has been implemented.")
        pf = self._mvr / 2.0
        egams = pf * xs

        # Check if `state` was passed in
        state = kwargs["state"] if "state" in kwargs else None
        if state is not None:
            if state in self._dnde_dispatch:
                return pf * self._dnde_dispatch[state](egams)
            raise ValueError(f"Invalid state: {state}")

        spectra = dict()
        for state, func in self._dnde_dispatch.items():
            spectra[state] = pf * func(egams)

        return spectra
