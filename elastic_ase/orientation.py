"""Crystallographic orientation: space group, whether a deformation table can represent the
crystal in its current orientation (``table_covers``), and rotations that align symmetry axes
with the Cartesian frame.
"""
from math import pi
import numpy as np

from elastic_ase import _spglib
from elastic_ase.auto import crystal_rotations, NOISE_TOL

SYMPREC = 1e-3
AXIS_TOL = 1e-5


def _spglib_cell(atoms):
    return (np.asarray(atoms.get_cell(), float), atoms.get_scaled_positions(), atoms.numbers)


def space_group(atoms, symprec=SYMPREC):
    """Return space group symbol and number string from spglib."""
    return _spglib.call("get_spacegroup", _spglib_cell(atoms), symprec=symprec) or "unknown"


def table_covers(spec, basis):
    """True if the named table's layout can represent every symmetry-allowed stiffness matrix."""
    n = len(spec.constants)
    table = np.array([spec.matrix(np.eye(n)[j]).ravel() for j in range(n)]).T
    for b in basis.values():
        target = np.asarray(b).ravel()
        coef, *_ = np.linalg.lstsq(table, target, rcond=None)
        if np.abs(table @ coef - target).max() > NOISE_TOL:
            return False
    return True


def _rotation(axis, angle):
    n = np.asarray(axis, float)
    n = n / np.linalg.norm(n)
    k = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _align(n, target):
    """Rotation matrix taking unit vector n onto target."""
    axis = np.cross(n, target)
    cos_angle = float(np.clip(n @ target, -1.0, 1.0))
    if np.linalg.norm(axis) < 1e-12:
        if cos_angle > 0:
            return np.eye(3)
        perp = np.cross(n, [1.0, 0.0, 0.0])
        if np.linalg.norm(perp) < 1e-6:
            perp = np.cross(n, [0.0, 1.0, 0.0])
        return _rotation(perp, pi)
    return _rotation(axis, np.arccos(cos_angle))


def symmetry_axes(atoms, order, symprec=SYMPREC):
    """Cartesian unit vectors of the crystal's n-fold rotation axes (n = order)."""
    try:
        rotations = crystal_rotations(atoms, symprec)
    except Exception:
        return []
    target_cos = np.cos(2 * pi / order)
    axes = []
    for r in rotations:
        if np.linalg.det(r) < 0:
            r = -r
        if abs((np.trace(r) - 1) / 2 - target_cos) > 1e-6:
            continue
        values, vectors = np.linalg.eig(r)
        axis = np.real(vectors[:, np.argmin(np.abs(values - 1))])
        axis /= np.linalg.norm(axis)
        if axis[np.argmax(np.abs(axis))] < 0:
            axis = -axis
        if not any(abs(abs(axis @ other) - 1) < AXIS_TOL for other in axes):
            axes.append(axis)
    return axes


def align_structure(atoms, order, target, symprec=SYMPREC):
    """
    Rotates the cell and Cartesian positions so that the nearest n-fold axis aligns with target.
    Preserves fractional coordinates, element order, and unit cell volume exactly.
    """
    axes = symmetry_axes(atoms, order, symprec)
    if not axes:
        return atoms, False

    target = np.asarray(target, float) / np.linalg.norm(target)
    # Check if already aligned
    if any(abs(abs(axis @ target) - 1) < AXIS_TOL for axis in axes):
        return atoms, False

    # Pick the axis closest to target
    best_axis = max(axes, key=lambda a: abs(a @ target))
    if best_axis @ target < 0:
        best_axis = -best_axis

    rot = _align(best_axis, target)
    rotated_atoms = atoms.copy()
    rotated_atoms.calc = atoms.calc
    # Cell transforms as cell @ rot^T
    rotated_atoms.set_cell(atoms.get_cell() @ rot.T, scale_atoms=False)
    # Positions transform as pos @ rot^T
    rotated_atoms.set_positions(atoms.get_positions() @ rot.T)
    return rotated_atoms, True
