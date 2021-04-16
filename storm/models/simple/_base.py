"""
Module containing the base class of all the `SimpleRhNeutrino` models.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

import numpy as np

from storm.constants import LEPTON_MASSES

StateType = Tuple[str, ...]


class SimpleRhNeutrinoBase(ABC):
    """
    Base class for all right-handed neutrino models where the right-handed
    neutrino mixes with a single active neutrino.
    """

    def __init__(self, mvr: float, theta: float, lep: str):
        self._mvr: float = mvr
        self._theta: float = theta
        self._lep: str = lep

        if lep not in ["e", "mu", "tau"]:
            raise ValueError(
                f"Invalid value for 'lep': {lep}. Use 'e', 'mu' or 'tau'")
        if lep == "e":
            self._genl = 0
        elif lep == "mu":
            self._genl = 1
        elif lep == "tau":
            self._genl = 2
        self._ml = LEPTON_MASSES[self._genl]

    @property
    def mvr(self) -> float:
        """
        The right-handed neutrino mass in GeV.
        """
        return self._mvr

    @mvr.setter
    def mvr(self, mvr: float) -> None:
        self._mvr = mvr

    @property
    def theta(self) -> float:
        """
        The mixing angle between the right-handed neutrino and the active
        neutrino.
        """
        return self._theta

    @theta.setter
    def theta(self, theta: float) -> None:
        self._theta = theta

    @property
    def lep(self) -> str:
        """
        String representing the flavor of neutrino that the right-handed
        neutrino mixes with.
        """
        return self._lep

    @lep.setter
    def lep(self, lep: str) -> None:
        self._lep = lep

    @property
    @abstractmethod
    def decay_final_states(self) -> List[StateType]:
        """
        Returns a list of all availible final states.
        """

    @abstractmethod
    def partial_width(self, state: StateType, **kwargs) -> float:
        """
        Compute the partial widths from the decay of a right-handed neutrino
        into a specified final state.
        """

    @abstractmethod
    def partial_widths(self, **kwargs) -> Dict[StateType, float]:
        """
        Compute the partial widths from the decay of a right-handed neutrino
        into all possible final states.
        """

    def branching_fractions(self, **kwargs) -> Dict[StateType, float]:
        """
        Compute the branching fractions from the decay of a right-handed
        neutrino into all possible final states.
        """
        widths = self.partial_widths()

        bfs = dict()
        for key in widths.keys():
            bfs[key] = widths[key] / widths[("total",)]

        return bfs

    @abstractmethod
    def dndx_single_state(
            self,
            x: np.ndarray,
            product: int,
            state: StateType,
            **kwargs
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the spectrum of a specified product from the decay of a
        right-handed neutrino into a specified final state.

        Parameters
        ----------
        scaled_energies: np.ndarray
            Array of the scaled product energies where the spectrum should be
            computed.
        product: int
            PDG code of the product to compute spectrum for.
        state: StateType
            Final state to compute spectrum for.

        Returns
        -------
        xs: np.ndarray
            Array of the x=2E/mvr values.
        dndx: np.ndarray
            Array of the spectrum.
        """

    @abstractmethod
    def dndx(
            self,
            x: np.ndarray,
            product: int,
            **kwargs
    ) -> Tuple[np.ndarray, Dict[StateType, np.ndarray]]:
        """
        Compute the spectrum of a specified product from the decay of a
        right-handed neutrino into all possible final states.

        Parameters
        ----------
        scaled_energies: np.ndarray
            Array of the scaled product energies where the spectrum should be
            computed.
        product: int
            PDG code of the product to compute spectrum for.

        Returns
        -------
        xs: np.ndarray
            Array of the x=2E/mvr values.
        spec_dict: Dict[StateType, np.ndarray]
            Dictionary containing spectra from all final states.
        """
