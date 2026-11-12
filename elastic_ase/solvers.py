import numpy as np

def check_mechanical_stability(stiffness, symmetry="cubic"):
    """
    Checks the Born mechanical stability criteria for the given stiffness matrix.
    """
    sym = symmetry.lower().strip()
    C = stiffness

    if sym == "cubic":
        C11, C12, C44 = C[0, 0], C[0, 1], C[3, 3]
        stable = (C11 - C12 > 0) and (C11 + 2.0 * C12 > 0) and (C44 > 0)
        criteria = "C11 - C12 > 0 and C11 + 2*C12 > 0 and C44 > 0"
    elif sym in ["hexagonal", "hex"]:
        C11, C12, C13, C33, C44 = C[0, 0], C[0, 1], C[0, 2], C[2, 2], C[3, 3]
        a = 2.0 * (C13 ** 2)
        b = C33 * (C11 + C12)
        stable = (C11 > abs(C12)) and (C44 > 0) and (C33 > 0) and ((C11 - C12) > 0) and (a < b)
        criteria = "C11 > |C12| and C44 > 0 and C33 > 0 and C11 - C12 > 0 and 2*C13^2 < C33*(C11 + C12)"
    else:
        # General positive-definite condition: all eigenvalues must be strictly positive
        C_sym = 0.5 * (C + C.T)
        eigenvalues = np.linalg.eigvalsh(C_sym)
        stable = bool(np.all(eigenvalues > 0))
        criteria = "All eigenvalues of stiffness matrix C > 0 (Positive Definite)"

    return stable, criteria

def solve_stiffness(symmetry, coeffs):
    """
    Assembles the 6x6 stiffness matrix (in GPa) from the fitted deformation coefficients.
    Formulas correspond to the strain-energy methodology in elastic_vasp.
    """
    sym = symmetry.lower().strip()

    if sym in ["cubic"]:
        c11_c12_I = coeffs["C11_C12_I"]
        c11_c12_II = coeffs["C11_C12_II"]
        c44 = coeffs["C44"]

        bulk = c11_c12_II * 2.0 / 9.0
        shear = c11_c12_I * 0.5
        C11 = (3.0 * bulk + 4.0 * shear) / 3.0
        C12 = (3.0 * bulk - 2.0 * shear) / 3.0
        C44 = 2.0 * c44 / 3.0

        C = np.array([
            [C11, C12, C12, 0.0, 0.0, 0.0],
            [C12, C11, C12, 0.0, 0.0, 0.0],
            [C12, C12, C11, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, C44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, C44, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, C44]
        ])
        return C

    elif sym in ["hexagonal", "hex"]:
        c11_c12_I = coeffs["C11_C12_I"]
        c11_c12_II = coeffs["C11_C12_II"]
        c11_c33_c13_I = coeffs["C11_C33_C13_I"]
        c44 = coeffs["C44"]
        c11_c12_c13_c33 = coeffs["C11_C12_C13_C33"]

        C11 = (c11_c12_I + c11_c12_II) / 2.0
        C12 = (c11_c12_I - c11_c12_II) / 2.0
        C44 = c44
        eq1 = 2.0 * c11_c33_c13_I - C11
        eq2 = 2.0 * (c11_c12_c13_c33 - c11_c12_I)
        C13 = (eq2 - eq1) / 6.0
        C33 = eq1 + 2.0 * C13
        C66 = (C11 - C12) / 2.0

        C = np.array([
            [C11, C12, C13, 0.0, 0.0, 0.0],
            [C12, C11, C13, 0.0, 0.0, 0.0],
            [C13, C13, C33, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, C44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, C44, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, C66]
        ])
        return C

    elif sym in ["orthogonal", "orthorhombic", "tetragonal"]:
        A = coeffs["C11_C22_C12_I"]
        B = coeffs["C11_C22_C12_II"]
        C_c = coeffs["C22_C33_C23_I"]
        D = coeffs["C22_C33_C23_II"]
        E = coeffs["C11_C33_C13_I"]
        F = coeffs["C11_C33_C13_II"]
        G = coeffs["C44"]
        H = coeffs["C55"]
        I = coeffs["C66"]

        AA = (A + B) / 2.0
        BB = (B - A) / 2.0
        CC = (D + C_c) / 2.0
        DD = (D - C_c) / 2.0
        EE = (E + F) / 2.0
        FF = (F - E) / 2.0

        C11 = AA - CC + EE
        C22 = AA + CC - EE
        C33 = -AA + CC + EE
        C12 = BB
        C23 = DD
        C13 = FF
        C44 = G * 2.0
        C55 = H * 2.0
        C66 = I * 2.0

        C = np.array([
            [C11, C12, C13, 0.0, 0.0, 0.0],
            [C12, C22, C23, 0.0, 0.0, 0.0],
            [C13, C23, C33, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, C44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, C55, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, C66]
        ])
        return C

    elif sym in ["trigonal_1", "trigonal1"]:
        A = coeffs["C11_C12_I"]
        B = coeffs["C11_C12_II"]
        C_c = coeffs["C11_C33_C13_I"]
        D = coeffs["C11_C33_C13_II"]
        E = coeffs["C11_C44_C14_I"]
        F = coeffs["C11_C44_C14_II"]

        C11 = (A + B) / 2.0
        C12 = (A - B) / 2.0
        C33 = (C_c + D) - C11
        C13 = (C_c - D) / 2.0
        C44 = (E + F) - C11
        C14 = (E - F) / 2.0
        C66 = (C11 - C12) / 2.0

        C = np.array([
            [C11, C12, C13, C14, 0.0, 0.0],
            [C12, C11, C13, -C14, 0.0, 0.0],
            [C13, C13, C33, 0.0, 0.0, 0.0],
            [C14, -C14, 0.0, C44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, C44, C14],
            [0.0, 0.0, 0.0, 0.0, C14, C66]
        ])
        return C

    elif sym in ["trigonal_2", "trigonal2"]:
        A = coeffs["C11_C12_I"]
        B = coeffs["C11_C12_II"]
        C_c = coeffs["C11_C33_C13_I"]
        D = coeffs["C11_C33_C13_II"]
        E = coeffs["C11_C44_C14_I"]
        F = coeffs["C11_C44_C14_II"]
        G = coeffs["C11_C44_C15"]

        C11 = (A + B) / 2.0
        C12 = (A - B) / 2.0
        C33 = (C_c + D) - C11
        C13 = (C_c - D) / 2.0
        C44 = (E + F) - C11
        C14 = (E - F) / 2.0
        C15 = (C11 + C44) / 2.0 - G
        C66 = (C11 - C12) / 2.0

        C = np.array([
            [C11, C12, C13, C14, C15, 0.0],
            [C12, C11, C13, -C14, -C15, 0.0],
            [C13, C13, C33, 0.0, 0.0, 0.0],
            [C14, -C14, 0.0, C44, 0.0, -C15],
            [C15, -C15, 0.0, 0.0, C44, C14],
            [0.0, 0.0, 0.0, -C15, C14, C66]
        ])
        return C

    elif sym in ["monoclinic"]:
        A = coeffs["C11_C22_C12_I"]
        B = coeffs["C11_C22_C12_II"]
        C_c = coeffs["C22_C33_C23_I"]
        D = coeffs["C22_C33_C23_II"]
        E = coeffs["C11_C33_C13_I"]
        F = coeffs["C11_C33_C13_II"]
        G = coeffs["C44"]
        H = coeffs["C11_C55_C15_I"]
        I = coeffs["C11_C55_C15_II"]
        J = coeffs["C22_C55_C25"]
        K = coeffs["C33_C55_C35"]
        L = coeffs["C44_C66_C46_I"]
        M = coeffs["C44_C66_C46_II"]

        AA = (A + B) / 2.0
        BB = (B - A) / 2.0
        CC = (D + C_c) / 2.0
        DD = (D - C_c) / 2.0
        EE = (E + F) / 2.0
        FF = (F - E) / 2.0

        C11 = AA - CC + EE
        C22 = AA + CC - EE
        C33 = -AA + CC + EE
        C12 = BB
        C23 = DD
        C13 = FF
        C44 = G * 2.0
        C55 = (H + I) - C11
        C15 = (H - I) / 2.0
        C25 = J - (C22 / 2.0) - (C55 / 2.0)
        C35 = K - (C33 / 2.0) - (C55 / 2.0)
        C66 = (L + M) - C44
        C46 = (L - M) / 2.0

        C = np.array([
            [C11, C12, C13, 0.0, C15, 0.0],
            [C12, C22, C23, 0.0, C25, 0.0],
            [C13, C23, C33, 0.0, C35, 0.0],
            [0.0, 0.0, 0.0, C44, 0.0, C46],
            [C15, C25, C35, 0.0, C55, 0.0],
            [0.0, 0.0, 0.0, C46, 0.0, C66]
        ])
        return C

    elif sym in ["triclinic"]:
        C11 = 2.0 * coeffs["C11"]
        C22 = 2.0 * coeffs["C22"]
        C33 = 2.0 * coeffs["C33"]
        C12 = coeffs["C11_C22_C12"] - (C11 + C22) / 2.0
        C13 = coeffs["C11_C33_C13"] - (C11 + C33) / 2.0
        C23 = coeffs["C22_C33_C23"] - (C22 + C33) / 2.0
        C44 = 2.0 * coeffs["C44"]
        C55 = 2.0 * coeffs["C55"]
        C66 = 2.0 * coeffs["C66"]
        C14 = coeffs["C11_C44_C14"] - (C11 + C44) / 2.0
        C15 = coeffs["C11_C55_C15"] - (C11 + C55) / 2.0
        C16 = coeffs["C11_C66_C16"] - (C11 + C66) / 2.0
        C24 = coeffs["C22_C44_C24"] - (C22 + C44) / 2.0
        C25 = coeffs["C22_C55_C25"] - (C22 + C55) / 2.0
        C26 = coeffs["C22_C66_C26"] - (C22 + C66) / 2.0
        C34 = coeffs["C33_C44_C34"] - (C33 + C44) / 2.0
        C35 = coeffs["C33_C55_C35"] - (C33 + C55) / 2.0
        C36 = coeffs["C33_C66_C36"] - (C33 + C66) / 2.0
        C45 = coeffs["C44_C55_C45"] - (C44 + C55) / 2.0
        C46 = coeffs["C44_C66_C46"] - (C44 + C66) / 2.0
        C56 = coeffs["C55_C66_C56"] - (C55 + C66) / 2.0

        C = np.array([
            [C11, C12, C13, C14, C15, C16],
            [C12, C22, C23, C24, C25, C26],
            [C13, C23, C33, C34, C35, C36],
            [C14, C24, C34, C44, C45, C46],
            [C15, C25, C35, C45, C55, C56],
            [C16, C26, C36, C46, C56, C66]
        ])
        return C

    else:
        raise ValueError(f"Unknown symmetry '{symmetry}'")
