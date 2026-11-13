"""Strain helpers. The deformation tables themselves live in :mod:`elastic_ase.systems`."""
import numpy as np

from elastic_ase.systems import get_system


def get_deformation_specs(symmetry):
    """[(folder, deformation)] for a crystal system, where deformation(delta) returns the Voigt
    strain (e1, ..., e6) of that deformation (kept for backward compatibility; the same data is
    in ``elastic_ase.systems.SYSTEMS[symmetry].deformations``)."""
    return [(d.folder, d.at) for d in get_system(symmetry).deformations]


def build_strain_matrix(d1, d2, d3, d4, d5, d6):
    """F = I + epsilon for a Voigt strain with engineering shear components (e4 = 2 eps_yz, ...)."""
    return np.array([
        [1.0 + d1, d6 / 2.0, d5 / 2.0],
        [d6 / 2.0, 1.0 + d2, d4 / 2.0],
        [d5 / 2.0, d4 / 2.0, 1.0 + d3],
    ])
