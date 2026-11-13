"""Crystal-system definitions and design-matrix solvers.

Each crystal system defines:
- its independent elastic constants,
- how they fill the 6x6 Voigt stiffness matrix, and
- the deformations: folder name, linear strain direction e (Voigt), and optional second-order correction.

For a deformation of size s, the strain energy density is E/V0 = a * s^2 where
a = 1/2 e^T C e (to leading order). The design matrix A has entries:
A_kj = 1/2 e_k^T C_j e_k
where C_j is the unit basis matrix for constant j. The independent constants
p are solved via linear least squares: A p = a.
"""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np


def _square_ratio(value, denominator):
    if abs(denominator) < 1e-15:
        raise ValueError(f"Invalid strain value {value} produces division by zero")
    return (value ** 2) / denominator


@dataclass(frozen=True)
class Deformation:
    folder: str
    strain: tuple
    # (component index, k): sets that component to delta^2 / (k - delta^2) for volume preservation
    correction: tuple = None

    def at(self, delta):
        """Voigt strain (e1..e6) applied for deformation size delta."""
        values = [0.0 if c == 0 else c * delta for c in self.strain]
        if self.correction is not None:
            index, k = self.correction
            values[index] = _square_ratio(delta, k - delta ** 2)
        return tuple(values)

    def describe(self):
        """Plain-text description of the strain."""
        text = "(" + ", ".join(f"{c:g}" for c in self.strain) + ") * delta"
        if self.correction is not None:
            index, k = self.correction
            text += f", e{index + 1} = delta^2/({k:g}-delta^2)"
        return text


@dataclass(frozen=True)
class CrystalSystem:
    name: str
    constants: tuple
    layout: object  # callable: namespace of constants -> 6x6 nested list
    deformations: tuple
    stability_text: str = "stiffness matrix is positive definite (all eigenvalues positive)"

    def matrix(self, values):
        """6x6 stiffness matrix from a sequence of independent constants."""
        ns = SimpleNamespace(**dict(zip(self.constants, values)))
        return np.array(self.layout(ns), dtype=float)

    def design_matrix(self):
        """Design matrix A where a = A p."""
        basis = np.eye(len(self.constants))
        unit_matrices = [self.matrix(row) for row in basis]
        design = np.empty((len(self.deformations), len(self.constants)))
        for k, deformation in enumerate(self.deformations):
            e = np.array(deformation.strain, dtype=float)
            for j, cj in enumerate(unit_matrices):
                design[k, j] = 0.5 * e @ cj @ e
        return design

    def energy_expression(self, k, design=None):
        """Human-readable a_k = 1/2 e^T C e for deformation k."""
        if design is None:
            design = self.design_matrix()
        terms = []
        for coeff, name in zip(design[k], self.constants):
            if abs(coeff) < 1e-12:
                continue
            sign = "-" if coeff < 0 else "+"
            mag = abs(coeff)
            body = name if abs(mag - 1) < 1e-12 else f"{mag:g}*{name}"
            terms.append((sign, body))
        if not terms:
            return "0"
        first_sign, first = terms[0]
        text = ("-" if first_sign == "-" else "") + first
        for sign, body in terms[1:]:
            text += f" {sign} {body}"
        return text

    def solve(self, coefficients):
        """Independent constants and 6x6 matrix from fitted coefficients."""
        design = self.design_matrix()
        rank = np.linalg.matrix_rank(design)
        if rank < len(self.constants):
            raise ValueError(
                f"{self.name}: deformations determine only {rank} of {len(self.constants)} constants"
            )
        rhs = np.array([coefficients[d.folder] for d in self.deformations], dtype=float)
        values, *_ = np.linalg.lstsq(design, rhs, rcond=None)
        return dict(zip(self.constants, values)), self.matrix(values)


def _d(folder, strain, correction=None):
    return Deformation(folder, tuple(strain), correction)


# --- Stiffness Layouts ---

def _cubic(c):
    return [[c.C11, c.C12, c.C12, 0, 0, 0],
            [c.C12, c.C11, c.C12, 0, 0, 0],
            [c.C12, c.C12, c.C11, 0, 0, 0],
            [0, 0, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C44, 0],
            [0, 0, 0, 0, 0, c.C44]]


def _hexagonal(c):
    c66 = (c.C11 - c.C12) / 2
    return [[c.C11, c.C12, c.C13, 0, 0, 0],
            [c.C12, c.C11, c.C13, 0, 0, 0],
            [c.C13, c.C13, c.C33, 0, 0, 0],
            [0, 0, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C44, 0],
            [0, 0, 0, 0, 0, c66]]


def _tetragonal(c, c16=0):
    return [[c.C11, c.C12, c.C13, 0, 0, c16],
            [c.C12, c.C11, c.C13, 0, 0, -c16],
            [c.C13, c.C13, c.C33, 0, 0, 0],
            [0, 0, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C44, 0],
            [c16, -c16, 0, 0, 0, c.C66]]


def _tetragonal_1(c):
    return _tetragonal(c)


def _tetragonal_2(c):
    return _tetragonal(c, c.C16)


def _orthorhombic(c):
    return [[c.C11, c.C12, c.C13, 0, 0, 0],
            [c.C12, c.C22, c.C23, 0, 0, 0],
            [c.C13, c.C23, c.C33, 0, 0, 0],
            [0, 0, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C55, 0],
            [0, 0, 0, 0, 0, c.C66]]


def _trigonal_1(c):
    c66 = (c.C11 - c.C12) / 2
    return [[c.C11, c.C12, c.C13, c.C14, 0, 0],
            [c.C12, c.C11, c.C13, -c.C14, 0, 0],
            [c.C13, c.C13, c.C33, 0, 0, 0],
            [c.C14, -c.C14, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C44, c.C14],
            [0, 0, 0, 0, c.C14, c66]]


def _trigonal_2(c):
    c66 = (c.C11 - c.C12) / 2
    return [[c.C11, c.C12, c.C13, c.C14, c.C15, 0],
            [c.C12, c.C11, c.C13, -c.C14, -c.C15, 0],
            [c.C13, c.C13, c.C33, 0, 0, 0],
            [c.C14, -c.C14, 0, c.C44, 0, -c.C15],
            [c.C15, -c.C15, 0, 0, c.C44, c.C14],
            [0, 0, 0, -c.C15, c.C14, c66]]


def _monoclinic(c):
    return [[c.C11, c.C12, c.C13, 0, c.C15, 0],
            [c.C12, c.C22, c.C23, 0, c.C25, 0],
            [c.C13, c.C23, c.C33, 0, c.C35, 0],
            [0, 0, 0, c.C44, 0, c.C46],
            [c.C15, c.C25, c.C35, 0, c.C55, 0],
            [0, 0, 0, c.C46, 0, c.C66]]


_TRICLINIC_NAMES = tuple(f"C{i}{j}" for i in range(1, 7) for j in range(i, 7))


def _triclinic(c):
    values = vars(c)
    return [[values[f"C{min(i, j)}{max(i, j)}"] for j in range(1, 7)] for i in range(1, 7)]


_ORTHORHOMBIC_NORMAL = (
    _d("C11_C22_C12_I", (1, -1, 0, 0, 0, 0), (2, 1)),
    _d("C11_C22_C12_II", (1, 1, 0, 0, 0, 0)),
    _d("C22_C33_C23_I", (0, 1, -1, 0, 0, 0), (0, 1)),
    _d("C22_C33_C23_II", (0, 1, 1, 0, 0, 0)),
    _d("C11_C33_C13_I", (1, 0, -1, 0, 0, 0), (1, 1)),
    _d("C11_C33_C13_II", (1, 0, 1, 0, 0, 0)),
    _d("C44", (0, 0, 0, 1, 0, 0), (0, 4)),
)

_TETRAGONAL_COMMON = (
    _d("C11_C12_I", (1, 1, 0, 0, 0, 0)),
    _d("C11_C12_II", (1, -1, 0, 0, 0, 0), (2, 1)),
    _d("C11_C33_C13_I", (1, 0, -1, 0, 0, 0), (1, 1)),
    _d("C11_C12_C13_C33", (1, 1, 1, 0, 0, 0)),
    _d("C44", (0, 0, 0, 1, 1, 0)),
    _d("C66", (0, 0, 0, 0, 0, 1), (2, 4)),
)

_TRIGONAL_COMMON = (
    _d("C11_C12_I", (1, 1, 0, 0, 0, 0)),
    _d("C11_C12_II", (1, -1, 0, 0, 0, 0), (2, 1)),
    _d("C11_C33_C13_I", (1, 0, 1, 0, 0, 0)),
    _d("C11_C33_C13_II", (1, 0, -1, 0, 0, 0), (1, 1)),
    _d("C11_C44_C14_I", (1, 0, 0, 1, 0, 0)),
    _d("C11_C44_C14_II", (1, 0, 0, -1, 0, 0)),
)

SYSTEMS = {
    "cubic": CrystalSystem(
        "cubic", ("C11", "C12", "C44"), _cubic,
        (
            _d("C11_C12_I", (1, -1, 0, 0, 0, 0)),
            _d("C11_C12_II", (1, 1, 1, 0, 0, 0)),
            _d("C44", (0, 0, 0, 1, 1, 1)),
        ),
        "C11-C12 > 0 ; C11+2C12 > 0 ; C44 > 0 (equivalent to positive definiteness)",
    ),
    "orthogonal": CrystalSystem(
        "orthogonal", ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66"),
        _orthorhombic,
        _ORTHORHOMBIC_NORMAL + (
            _d("C55", (0, 0, 0, 0, 1, 0), (1, 4)),
            _d("C66", (0, 0, 0, 0, 0, 1), (2, 4)),
        ),
    ),
    "hexagonal": CrystalSystem(
        "hexagonal", ("C11", "C12", "C13", "C33", "C44"), _hexagonal,
        (
            _d("C11_C12_I", (1, 1, 0, 0, 0, 0)),
            _d("C11_C12_II", (1, -1, 0, 0, 0, 0), (2, 1)),
            _d("C11_C33_C13_I", (1, 0, -1, 0, 0, 0), (1, 1)),
            _d("C44", (0, 0, 0, 1, 1, 0)),
            _d("C11_C12_C13_C33", (1, 1, 1, 0, 0, 0)),
        ),
        "C11 > |C12| ; 2C13^2 < C33(C11+C12) ; C44 > 0 ; C66 > 0",
    ),
    "tetragonal_1": CrystalSystem(
        "tetragonal_1", ("C11", "C12", "C13", "C33", "C44", "C66"), _tetragonal_1,
        _TETRAGONAL_COMMON,
    ),
    "tetragonal_2": CrystalSystem(
        "tetragonal_2", ("C11", "C12", "C13", "C33", "C44", "C66", "C16"), _tetragonal_2,
        _TETRAGONAL_COMMON + (_d("C11_C66_C16", (1, 0, 0, 0, 0, 1)),),
    ),
    "trigonal_1": CrystalSystem(
        "trigonal_1", ("C11", "C12", "C13", "C33", "C44", "C14"), _trigonal_1,
        _TRIGONAL_COMMON,
    ),
    "trigonal_2": CrystalSystem(
        "trigonal_2", ("C11", "C12", "C13", "C33", "C44", "C14", "C15"), _trigonal_2,
        _TRIGONAL_COMMON + (_d("C11_C44_C15", (1, 0, 0, 0, -1, 0)),),
    ),
    "monoclinic": CrystalSystem(
        "monoclinic",
        ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66",
         "C15", "C25", "C35", "C46"),
        _monoclinic,
        _ORTHORHOMBIC_NORMAL + (
            _d("C11_C55_C15_I", (1, 0, 0, 0, 1, 0)),
            _d("C11_C55_C15_II", (1, 0, 0, 0, -1, 0)),
            _d("C22_C55_C25", (0, 1, 0, 0, 1, 0)),
            _d("C33_C55_C35", (0, 0, 1, 0, 1, 0)),
            _d("C44_C66_C46_I", (0, 0, 0, 1, 0, 1)),
            _d("C44_C66_C46_II", (0, 0, 0, 1, 0, -1)),
        ),
    ),
    "triclinic": CrystalSystem(
        "triclinic", _TRICLINIC_NAMES, _triclinic,
        (
            _d("C11", (1, 0, 0, 0, 0, 0)),
            _d("C11_C22_C12", (1, 1, 0, 0, 0, 0)),
            _d("C22", (0, 1, 0, 0, 0, 0)),
            _d("C22_C33_C23", (0, 1, 1, 0, 0, 0)),
            _d("C33", (0, 0, 1, 0, 0, 0)),
            _d("C11_C33_C13", (1, 0, 1, 0, 0, 0)),
            _d("C44", (0, 0, 0, 1, 0, 0)),
            _d("C55", (0, 0, 0, 0, 1, 0)),
            _d("C66", (0, 0, 0, 0, 0, 1)),
            _d("C11_C44_C14", (1, 0, 0, 1, 0, 0)),
            _d("C11_C55_C15", (1, 0, 0, 0, 1, 0)),
            _d("C11_C66_C16", (1, 0, 0, 0, 0, 1)),
            _d("C22_C44_C24", (0, 1, 0, 1, 0, 0)),
            _d("C22_C55_C25", (0, 1, 0, 0, 1, 0)),
            _d("C22_C66_C26", (0, 1, 0, 0, 0, 1)),
            _d("C33_C44_C34", (0, 0, 1, 1, 0, 0)),
            _d("C33_C55_C35", (0, 0, 1, 0, 1, 0)),
            _d("C33_C66_C36", (0, 0, 1, 0, 0, 1)),
            _d("C44_C55_C45", (0, 0, 0, 1, 1, 0)),
            _d("C44_C66_C46", (0, 0, 0, 1, 0, 1)),
            _d("C55_C66_C56", (0, 0, 0, 0, 1, 1)),
        ),
    ),
}

# Aliases
SYSTEMS["orthorhombic"] = SYSTEMS["orthogonal"]

_ALIASES = {"hex": "hexagonal", "tetragonal": "tetragonal_1", "tetragonal1": "tetragonal_1",
            "tetragonal2": "tetragonal_2", "trigonal1": "trigonal_1", "trigonal2": "trigonal_2"}


def get_system(name):
    """The CrystalSystem for a table name (or a common alias such as "orthorhombic")."""
    key = str(name).lower().strip()
    key = _ALIASES.get(key, key)
    if key not in SYSTEMS:
        raise ValueError(f"Unsupported symmetry: {name!r}. Choose from {sorted(set(SYSTEMS) - {'orthorhombic'})}.")
    return SYSTEMS[key]
