"""
Module containing the version of the `SimpleRhNeutrino` model valid for
right-handed neutrino masses less than 1 GeV.
"""

from typing import Callable, Dict, List, Tuple

import numpy as np
from hazma.rh_neutrino import RHNeutrino as _HazmaRhNeutrino  # type: ignore

from storm.models.simple._base import SimpleRhNeutrinoBase, StateType


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
            Callable[..., np.ndarray]
        ] = dict()
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

    def dndx_single_state(
            self,
            x: np.ndarray,
            product: int,
            state: StateType,
            **kwargs
    ) -> Tuple[np.ndarray, np.ndarray]:
        if not product == "photon" or not product == 22:
            raise NotImplementedError(
                "Only the photon spectrum has been implemented.")
        pf = self._mvr / 2.0
        egams = pf * x

        # Check if `state` was passed in
        if state in self._dnde_dispatch:
            return x, pf * self._dnde_dispatch[state](egams)
        raise ValueError(f"Invalid state: {state}")

    def dndx(
            self,
            x: np.ndarray,
            product: int,
            **kwargs
    ) -> Tuple[np.ndarray, Dict[StateType, np.ndarray]]:
        if not product == "photon" or not product == 22:
            raise NotImplementedError(
                "Only the photon spectrum has been implemented.")
        pf = self._mvr / 2.0
        egams = pf * x

        spectra = dict()
        for state, func in self._dnde_dispatch.items():
            spectra[state] = pf * func(egams)

        return x, spectra
