"""Orientation-aware automatic mode: independent constants and deformations from the crystal.

Instead of assuming a standard orientation (c along z, unique axis along y, etc.),
the crystal's symmetry operations are derived via spglib from the primitive cell
in the user's own Cartesian frame:

1. The stiffness matrices invariant under every point group operation form a linear
   subspace. Its reduced row-echelon basis over C11..C66 identifies the
   independent constants in this exact orientation.
2. A minimal, well-conditioned set of deformations is selected greedily from simple
   candidates, maximizing the smallest singular value (sigma_min).
3. The resulting stiffness tensor can be reported both in the user's frame and
   transformed to the textbook standard crystallographic setting.
"""
from itertools import combinations
import numpy as np

from elastic_ase._spglib import call as spg_call
from elastic_ase.systems import CrystalSystem, Deformation

SYMPREC = 1e-3
_TOL = 1e-8
NOISE_TOL = 1e-6

# Voigt index -> Cartesian index pair
_VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
_UPPER = [(i, j) for i in range(6) for j in range(i, 6)]
_NAMES = [f"C{i + 1}{j + 1}" for i, j in _UPPER]


def bond_matrix(rotation):
    """6x6 Bond matrix M with C' = M C M^T for a Cartesian rotation."""
    r = np.asarray(rotation, float)
    m = np.zeros((6, 6))
    for a, (i, j) in enumerate(_VOIGT):
        for b, (k, l) in enumerate(_VOIGT):
            m[a, b] = r[i, k] * r[j, l] if b < 3 else r[i, k] * r[j, l] + r[i, l] * r[j, k]
    return m


def _to_matrix(vector):
    c = np.zeros((6, 6))
    for value, (i, j) in zip(vector, _UPPER):
        c[i, j] = c[j, i] = value
    return c


def _to_vector(matrix):
    return np.array([matrix[i, j] for i, j in _UPPER])


def primitive_cell(atoms, symprec=SYMPREC):
    """spglib cell tuple (lattice, positions, numbers) of primitive cell, same orientation."""
    cell = (np.asarray(atoms.get_cell(), float), atoms.get_scaled_positions(), atoms.numbers)
    primitive = spg_call("standardize_cell", cell, to_primitive=True, no_idealize=True, symprec=symprec)
    return primitive if primitive is not None else cell


def crystal_rotations(atoms, symprec=SYMPREC):
    """Cartesian rotation matrices of the crystal's point group (from spglib)."""
    cell = primitive_cell(atoms, symprec)
    data = spg_call("get_symmetry", cell, symprec=symprec)
    if not data:
        raise ValueError(f"spglib could not determine symmetry (symprec={symprec})")
    lattice = np.asarray(cell[0])
    to_cart, from_cart = lattice.T, np.linalg.inv(lattice.T)
    rotations = []
    for w in data["rotations"]:
        r = to_cart @ w @ from_cart
        u, _, vt = np.linalg.svd(r)
        r = u @ vt
        if not any(np.allclose(r, other, atol=1e-6) for other in rotations):
            rotations.append(r)
    return rotations


def standard_setting(atoms, symprec=SYMPREC):
    """(rotation, space_group_number) taking this crystal to its standard setting."""
    dataset = spg_call("get_symmetry_dataset", primitive_cell(atoms, symprec), symprec=symprec)
    if dataset is None:
        raise ValueError(f"spglib could not determine symmetry (symprec={symprec})")
    get = (lambda key: getattr(dataset, key)) if hasattr(dataset, "number") else dataset.__getitem__
    rotation = np.array(get("std_rotation_matrix"), float)
    u, _, vt = np.linalg.svd(rotation)
    return u @ vt, int(get("number"))


def to_standard(c, rotation):
    """Stiffness matrix C expressed in the standard crystallographic frame."""
    m = bond_matrix(rotation)
    return m @ np.asarray(c, float) @ m.T


def invariant_basis(rotations):
    """{name: 6x6 matrix}: reduced row-echelon basis of symmetry-allowed stiffness matrices."""
    eye = np.eye(len(_UPPER))
    blocks = []
    for r in rotations:
        m = bond_matrix(r)
        t = np.array([_to_vector(m @ _to_matrix(e) @ m.T) for e in eye]).T  # 21x21
        blocks.append(t - eye)
    a = np.vstack(blocks) if blocks else np.zeros((1, len(_UPPER)))
    _, s, vt = np.linalg.svd(a)
    s = np.concatenate([s, np.zeros(vt.shape[0] - len(s))])
    null = vt[s < 1e-6 * max(1.0, s.max())]
    rows = _rref(null)
    basis = {}
    for row in rows:
        pivot = int(np.argmax(np.abs(row) > _TOL))
        row = np.where(np.abs(row) < 1e-10, 0.0, row)
        basis[_NAMES[pivot]] = _to_matrix(row)
    return basis


def _rref(m):
    m = np.array(m, float)
    rows, cols = m.shape
    lead = 0
    for c in range(cols):
        if lead >= rows:
            break
        p = lead + int(np.argmax(np.abs(m[lead:, c])))
        if abs(m[p, c]) < _TOL:
            continue
        m[[lead, p]] = m[[p, lead]]
        m[lead] /= m[lead, c]
        for r in range(rows):
            if r != lead:
                m[r] -= m[r, c] * m[lead]
        lead += 1
    return m[:lead]


def _candidates():
    """Candidate deformations in order of preference, with volume-conserving corrections."""
    out = []
    for i in range(6):
        e = [0] * 6
        e[i] = 1
        correction = (i - 3, 4) if i >= 3 else None
        out.append((e, correction))
    for i, j in combinations(range(6), 2):
        for sign in (1, -1):
            e = [0] * 6
            e[i], e[j] = 1, sign
            correction = None
            if sign < 0 and j < 3:
                correction = (3 - i - j, 1)
            out.append((e, correction))
    return out


def _folder_name(strain):
    parts = []
    for index, value in enumerate(strain):
        if value:
            label = f"e{index + 1}"
            parts.append(label if not parts else ("+" if value > 0 else "-") + label)
    return "".join(parts)


def _energy_row(strain, basis_matrices):
    e = np.asarray(strain, float)
    return np.array([0.5 * e @ b @ e for b in basis_matrices])


def _orthonormal(basis):
    vectors = np.array([np.asarray(m).ravel() for m in basis.values()]).T
    q, _ = np.linalg.qr(vectors)
    return [q[:, j].reshape(6, 6) for j in range(q.shape[1])]


def choose_deformations(basis):
    """Minimal set of deformations determining stiffness matrix with optimal condition number."""
    matrices = _orthonormal(basis)
    n = len(matrices)
    chosen, rows = [], []
    pool = _candidates()
    while len(rows) < n:
        best = None
        for strain, correction in pool:
            if any(strain == c[0] for c in chosen):
                continue
            e = np.asarray(strain, float)
            row = _energy_row(strain, matrices) / (e @ e)
            if np.linalg.norm(row) < _TOL:
                continue
            s = np.linalg.svd(np.vstack(rows + [row]), compute_uv=False)
            rank = int(np.sum(s > 1e-8))
            score = (rank, s[rank - 1] if rank else 0.0)
            if best is None or score[0] > best[0][0] or (
                    score[0] == best[0][0] and score[1] > best[0][1] * (1 + 1e-9)):
                best = (score, strain, correction, row)
        if best is None or best[0][0] <= len(rows):
            raise ValueError("Could not find deformations that determine all constants")
        chosen.append((best[1], best[2]))
        rows.append(best[3])
    return tuple(Deformation(_folder_name(s), tuple(s), c) for s, c in chosen)


def _layout(basis):
    names = list(basis)

    def layout(ns):
        return sum(getattr(ns, name) * basis[name] for name in names).tolist()
    return layout


def from_structure(atoms, symprec=SYMPREC):
    """Generates an orientation-adapted CrystalSystem for this crystal in its current Cartesian frame."""
    sg = spg_call("get_spacegroup", (np.asarray(atoms.get_cell(), float), atoms.get_scaled_positions(), atoms.numbers), symprec=symprec) or "unknown"
    basis = invariant_basis(crystal_rotations(atoms, symprec))
    spec = CrystalSystem(
        f"auto: {sg}",
        tuple(basis),
        _layout(basis),
        choose_deformations(basis)
    )
    return spec, basis


def standard_constants(c_standard, spacegroup_number):
    """Named table and textbook constants for a stiffness matrix in the standard setting."""
    from elastic_ase.symmetry import classify
    from elastic_ase.systems import SYSTEMS

    name = classify(spacegroup_number)[2]
    table = SYSTEMS[name]
    n = len(table.constants)
    design = np.array([table.matrix(np.eye(n)[j]).ravel() for j in range(n)]).T
    values, *_ = np.linalg.lstsq(design, np.asarray(c_standard).ravel(), rcond=None)
    residual = float(np.abs(design @ values - np.asarray(c_standard).ravel()).max())
    return name, dict(zip(table.constants, values)), residual
