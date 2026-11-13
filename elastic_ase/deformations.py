import numpy as np

def square_ratio(value, denominator):
    if abs(denominator) < 1e-15:
        raise ValueError(f"Strain value {value} produces division by zero")
    return (value ** 2) / denominator

def get_deformation_specs(symmetry):
    """
    Returns list of (folder_name, deformation_func) for the given crystal symmetry.
    Each deformation_func maps a strain scalar s -> (d1, d2, d3, d4, d5, d6) in Voigt notation.
    """
    sym = symmetry.lower().strip()

    if sym in ["cubic"]:
        return [
            ("C11_C12_I", lambda s: (s, -s, 0.0, 0.0, 0.0, 0.0)),
            ("C11_C12_II", lambda s: (s, s, s, 0.0, 0.0, 0.0)),
            ("C44", lambda s: (0.0, 0.0, 0.0, s, s, s)),
        ]
    elif sym in ["hexagonal", "hex"]:
        return [
            ("C11_C12_I", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
            ("C11_C12_II", lambda s: (s, -s, square_ratio(s, 1.0 - s**2), 0.0, 0.0, 0.0)),
            ("C11_C33_C13_I", lambda s: (s, square_ratio(s, 1.0 - s**2), -s, 0.0, 0.0, 0.0)),
            ("C44", lambda s: (0.0, 0.0, 0.0, s, s, 0.0)),
            ("C11_C12_C13_C33", lambda s: (s, s, s, 0.0, 0.0, 0.0)),
        ]
    elif sym in ["orthogonal", "orthorhombic"]:
        return [
            ("C11_C22_C12_I", lambda s: (s, -s, square_ratio(s, 1.0 - s**2), 0.0, 0.0, 0.0)),
            ("C11_C22_C12_II", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
            ("C22_C33_C23_I", lambda s: (square_ratio(s, 1.0 - s**2), s, -s, 0.0, 0.0, 0.0)),
            ("C22_C33_C23_II", lambda s: (0.0, s, s, 0.0, 0.0, 0.0)),
            ("C11_C33_C13_I", lambda s: (s, square_ratio(s, 1.0 - s**2), -s, 0.0, 0.0, 0.0)),
            ("C11_C33_C13_II", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
            ("C44", lambda s: (square_ratio(s, 4.0 - s**2), 0.0, 0.0, s, 0.0, 0.0)),
            ("C55", lambda s: (0.0, square_ratio(s, 4.0 - s**2), 0.0, 0.0, s, 0.0)),
            ("C66", lambda s: (0.0, 0.0, square_ratio(s, 4.0 - s**2), 0.0, 0.0, s)),
        ]
    elif sym in ["tetragonal_1", "tetragonal1", "tetragonal"]:
        return [
            ("C11_C12_I", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
            ("C11_C12_II", lambda s: (s, -s, square_ratio(s, 1.0 - s**2), 0.0, 0.0, 0.0)),
            ("C11_C33_C13_I", lambda s: (s, square_ratio(s, 1.0 - s**2), -s, 0.0, 0.0, 0.0)),
            ("C11_C12_C13_C33", lambda s: (s, s, s, 0.0, 0.0, 0.0)),
            ("C44", lambda s: (0.0, 0.0, 0.0, s, s, 0.0)),
            ("C66", lambda s: (0.0, 0.0, square_ratio(s, 4.0 - s**2), 0.0, 0.0, s)),
        ]
    elif sym in ["tetragonal_2", "tetragonal2"]:
        return [
            ("C11_C12_I", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
            ("C11_C12_II", lambda s: (s, -s, square_ratio(s, 1.0 - s**2), 0.0, 0.0, 0.0)),
            ("C11_C33_C13_I", lambda s: (s, square_ratio(s, 1.0 - s**2), -s, 0.0, 0.0, 0.0)),
            ("C11_C12_C13_C33", lambda s: (s, s, s, 0.0, 0.0, 0.0)),
            ("C44", lambda s: (0.0, 0.0, 0.0, s, s, 0.0)),
            ("C66", lambda s: (0.0, 0.0, square_ratio(s, 4.0 - s**2), 0.0, 0.0, s)),
            ("C11_C66_C16", lambda s: (s, 0.0, 0.0, 0.0, 0.0, s)),
        ]
    elif sym in ["trigonal_1", "trigonal1"]:
        return [
            ("C11_C12_I", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
            ("C11_C12_II", lambda s: (s, -s, square_ratio(s, 1.0 - s**2), 0.0, 0.0, 0.0)),
            ("C11_C33_C13_I", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
            ("C11_C33_C13_II", lambda s: (s, square_ratio(s, 1.0 - s**2), -s, 0.0, 0.0, 0.0)),
            ("C11_C44_C14_I", lambda s: (s, 0.0, 0.0, s, 0.0, 0.0)),
            ("C11_C44_C14_II", lambda s: (s, 0.0, 0.0, -s, 0.0, 0.0)),
        ]
    elif sym in ["trigonal_2", "trigonal2"]:
        return [
            ("C11_C12_I", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
            ("C11_C12_II", lambda s: (s, -s, square_ratio(s, 1.0 - s**2), 0.0, 0.0, 0.0)),
            ("C11_C33_C13_I", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
            ("C11_C33_C13_II", lambda s: (s, square_ratio(s, 1.0 - s**2), -s, 0.0, 0.0, 0.0)),
            ("C11_C44_C14_I", lambda s: (s, 0.0, 0.0, s, 0.0, 0.0)),
            ("C11_C44_C14_II", lambda s: (s, 0.0, 0.0, -s, 0.0, 0.0)),
            ("C11_C44_C15", lambda s: (s, 0.0, 0.0, 0.0, -s, 0.0)),
        ]
    elif sym in ["monoclinic"]:
        return [
            ("C11_C22_C12_I", lambda s: (s, -s, square_ratio(s, 1.0 - s**2), 0.0, 0.0, 0.0)),
            ("C11_C22_C12_II", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
            ("C22_C33_C23_I", lambda s: (square_ratio(s, 1.0 - s**2), s, -s, 0.0, 0.0, 0.0)),
            ("C22_C33_C23_II", lambda s: (0.0, s, s, 0.0, 0.0, 0.0)),
            ("C11_C33_C13_I", lambda s: (s, square_ratio(s, 1.0 - s**2), -s, 0.0, 0.0, 0.0)),
            ("C11_C33_C13_II", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
            ("C44", lambda s: (square_ratio(s, 4.0 - s**2), 0.0, 0.0, s, 0.0, 0.0)),
            ("C11_C55_C15_I", lambda s: (s, 0.0, 0.0, 0.0, s, 0.0)),
            ("C11_C55_C15_II", lambda s: (s, 0.0, 0.0, 0.0, -s, 0.0)),
            ("C22_C55_C25", lambda s: (0.0, s, 0.0, 0.0, s, 0.0)),
            ("C33_C55_C35", lambda s: (0.0, 0.0, s, 0.0, s, 0.0)),
            ("C44_C66_C46_I", lambda s: (0.0, 0.0, 0.0, s, 0.0, s)),
            ("C44_C66_C46_II", lambda s: (0.0, 0.0, 0.0, s, 0.0, -s)),
        ]
    elif sym in ["triclinic"]:
        return [
            ("C11", lambda s: (s, 0.0, 0.0, 0.0, 0.0, 0.0)),
            ("C11_C22_C12", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
            ("C22", lambda s: (0.0, s, 0.0, 0.0, 0.0, 0.0)),
            ("C22_C33_C23", lambda s: (0.0, s, s, 0.0, 0.0, 0.0)),
            ("C33", lambda s: (0.0, 0.0, s, 0.0, 0.0, 0.0)),
            ("C11_C33_C13", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
            ("C44", lambda s: (0.0, 0.0, 0.0, s, 0.0, 0.0)),
            ("C55", lambda s: (0.0, 0.0, 0.0, 0.0, s, 0.0)),
            ("C66", lambda s: (0.0, 0.0, 0.0, 0.0, 0.0, s)),
            ("C11_C44_C14", lambda s: (s, 0.0, 0.0, s, 0.0, 0.0)),
            ("C11_C55_C15", lambda s: (s, 0.0, 0.0, 0.0, s, 0.0)),
            ("C11_C66_C16", lambda s: (s, 0.0, 0.0, 0.0, 0.0, s)),
            ("C22_C44_C24", lambda s: (0.0, s, 0.0, s, 0.0, 0.0)),
            ("C22_C55_C25", lambda s: (0.0, s, 0.0, 0.0, s, 0.0)),
            ("C22_C66_C26", lambda s: (0.0, s, 0.0, 0.0, 0.0, s)),
            ("C33_C44_C34", lambda s: (0.0, 0.0, s, s, 0.0, 0.0)),
            ("C33_C55_C35", lambda s: (0.0, 0.0, s, 0.0, s, 0.0)),
            ("C33_C66_C36", lambda s: (0.0, 0.0, s, 0.0, 0.0, s)),
            ("C44_C55_C45", lambda s: (0.0, 0.0, 0.0, s, s, 0.0)),
            ("C44_C66_C46", lambda s: (0.0, 0.0, 0.0, s, 0.0, s)),
            ("C55_C66_C56", lambda s: (0.0, 0.0, 0.0, 0.0, s, s)),
        ]
    else:
        raise ValueError(f"Unsupported symmetry: '{symmetry}'. Choose from 'cubic', 'hexagonal', 'orthogonal', 'trigonal_1', 'trigonal_2', 'monoclinic', 'triclinic'.")

def build_strain_matrix(d1, d2, d3, d4, d5, d6):
    """
    Constructs the 3x3 strain matrix (I + epsilon) from Voigt strains.
    """
    return np.array([
        [1.0 + d1, d6 / 2.0, d5 / 2.0],
        [d6 / 2.0, 1.0 + d2, d4 / 2.0],
        [d5 / 2.0, d4 / 2.0, 1.0 + d3]
    ])
