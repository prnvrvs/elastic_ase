"""The results of one elastic-constant analysis."""
from collections.abc import Mapping
from dataclasses import dataclass, field, fields
from typing import Optional

import numpy as np


@dataclass
class ElasticResults(Mapping):
    """Everything computed from the energies. Moduli in GPa, energies in eV.

    Works like a read-only dictionary as well (``results["stiffness_matrix"]``,
    ``results.items()``), so code written for the earlier dictionary keeps working.
    """
    stiffness_matrix: np.ndarray          # in the frame chosen with frame=
    stiffness_cartesian: np.ndarray       # frame of the calculation (after any rotation)
    stiffness_standard: np.ndarray        # standard crystallographic setting
    standard_constants: dict              # independent constants of the space group, standard setting
    standard_rotation: Optional[np.ndarray]
    independent_constants: dict           # constants of the deformation table
    is_stable: bool                       # Born stability (all eigenvalues of C positive)
    stability_criteria: str
    polycrystalline_averages: dict        # Voigt / Reuss / Hill moduli and anisotropy indices
    directional_properties: Optional[dict]  # extremes from elastic-vis, None without it
    derived_properties: dict              # density, sound velocities, Debye temperature, ...
    residual_stress: dict                 # residual pressure of the reference cell from the fits
    coefficients: dict                    # fitted curvature a of every deformation (GPa)
    r2_scores: dict                       # R^2 of every quadratic fit
    energies: dict                        # every energy (eV) by run name
    fits: dict = field(default_factory=dict)  # strains, (E-E0)/V0 and fit polynomial per deformation

    # Mapping interface ------------------------------------------------------------------
    def __getitem__(self, key):
        if key in self._names():
            return getattr(self, key)
        raise KeyError(key)

    def __iter__(self):
        return iter(self._names())

    def __len__(self):
        return len(self._names())

    @classmethod
    def _names(cls):
        return [f.name for f in fields(cls)]

    def to_dict(self):
        """A plain dictionary (arrays kept as numpy arrays)."""
        return {name: getattr(self, name) for name in self._names()}
