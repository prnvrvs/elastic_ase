import numpy as np

from elastic_ase.systems import get_system

def check_mechanical_stability(stiffness, symmetry="cubic"):
    """
    Born mechanical stability: the stiffness matrix must be positive definite.

    The decision always uses the eigenvalues, which is exact in any orientation. For cubic and
    hexagonal crystals the equivalent textbook criteria are reported as the description; they
    hold only for a matrix in the standard setting (crystal axes along x, y, z).
    """
    sym = symmetry.lower().strip()
    C = np.asarray(stiffness, float)
    eigenvalues = np.linalg.eigvalsh(0.5 * (C + C.T))
    stable = bool(np.all(eigenvalues > 0))

    if sym == "cubic":
        criteria = "C11 - C12 > 0 and C11 + 2*C12 > 0 and C44 > 0"
    elif sym in ["hexagonal", "hex"]:
        criteria = "C11 > |C12| and C44 > 0 and C33 > 0 and C11 - C12 > 0 and 2*C13^2 < C33*(C11 + C12)"
    else:
        criteria = "All eigenvalues of stiffness matrix C > 0 (Positive Definite)"
    criteria += f" (smallest eigenvalue {eigenvalues.min():.2f} GPa)"

    return stable, criteria

def solve_stiffness(symmetry, coeffs):
    """6x6 stiffness matrix (GPa) from the fitted curvatures {deformation folder: a (GPa)}
    (kept for backward compatibility; same as ``SYSTEMS[symmetry].solve(coeffs)[1]``)."""
    return get_system(symmetry).solve(coeffs)[1]
