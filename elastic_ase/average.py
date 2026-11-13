import numpy as np

def calculate_polycrystalline_averages(C):
    """
    Compute Voigt, Reuss, and Hill polycrystalline elastic averages from 6x6 stiffness matrix C (in GPa).
    """
    C11 = C[0, 0]
    C22 = C[1, 1]
    C33 = C[2, 2]
    C44 = C[3, 3]
    C55 = C[4, 4]
    C66 = C[5, 5]
    C12 = C[0, 1]
    C13 = C[0, 2]
    C23 = C[1, 2]

    # Voigt averages
    K_voigt = ((C11 + C22 + C33) + 2.0 * (C12 + C13 + C23)) / 9.0
    G_voigt = ((C11 + C22 + C33) - (C12 + C13 + C23) + 3.0 * (C44 + C55 + C66)) / 15.0

    try:
        S = np.linalg.inv(C)
        S11, S22, S33 = S[0, 0], S[1, 1], S[2, 2]
        S44, S55, S66 = S[3, 3], S[4, 4], S[5, 5]
        S12, S13, S23 = S[0, 1], S[0, 2], S[1, 2]

        k_denom = (S11 + S22 + S33) + 2.0 * (S12 + S13 + S23)
        g_denom = 4.0 * (S11 + S22 + S33) - 4.0 * (S12 + S13 + S23) + 3.0 * (S44 + S55 + S66)

        if abs(k_denom) < 1e-15 or abs(g_denom) < 1e-15:
            K_reuss = K_voigt
            G_reuss = G_voigt
        else:
            K_reuss = 1.0 / k_denom
            G_reuss = 15.0 / g_denom
    except np.linalg.LinAlgError:
        K_reuss = K_voigt
        G_reuss = G_voigt

    K_VRH = 0.5 * (K_voigt + K_reuss)
    G_VRH = 0.5 * (G_voigt + G_reuss)

    denom_E = G_VRH + 3.0 * K_VRH
    E_VRH = (9.0 * K_VRH * G_VRH) / denom_E if abs(denom_E) > 1e-12 else 0.0

    denom_nu = 6.0 * K_VRH + 2.0 * G_VRH
    nu_VRH = (3.0 * K_VRH - 2.0 * G_VRH) / denom_nu if abs(denom_nu) > 1e-12 else 0.0

    pugh = K_VRH / G_VRH if abs(G_VRH) > 1e-12 else 0.0

    # Cauchy pressure and Zener anisotropy (for cubic/orthorhombic systems)
    cauchy = C12 - C44
    denom_A = C11 - C12
    anisotropy = (2.0 * C44) / denom_A if abs(denom_A) > 1e-12 else 1.0

    # Universal elastic anisotropy index A^U (Ranganathan & Ostoja-Starzewski, PRL 2008)
    # A^U = 5 * (G_voigt / G_reuss) + (K_voigt / K_reuss) - 6 >= 0 (0 for isotropic)
    if G_reuss > 1e-6 and K_reuss > 1e-6:
        universal_anisotropy = max(0.0, float(5.0 * (G_voigt / G_reuss) + (K_voigt / K_reuss) - 6.0))
    else:
        universal_anisotropy = 0.0

    # Chung-Buessem fractional anisotropy
    cb_G = (G_voigt - G_reuss) / (G_voigt + G_reuss) if (G_voigt + G_reuss) > 1e-6 else 0.0
    cb_K = (K_voigt - K_reuss) / (K_voigt + K_reuss) if (K_voigt + K_reuss) > 1e-6 else 0.0

    return {
        "K_voigt": float(K_voigt),
        "K_reuss": float(K_reuss),
        "K_VRH": float(K_VRH),
        "G_voigt": float(G_voigt),
        "G_reuss": float(G_reuss),
        "G_VRH": float(G_VRH),
        "E_VRH": float(E_VRH),
        "nu_VRH": float(nu_VRH),
        "Pugh_ratio": float(pugh),
        "Cauchy_pressure": float(cauchy),
        "Zener_anisotropy": float(anisotropy),
        "universal_anisotropy": float(universal_anisotropy),
        "chung_buessem_G": float(cb_G),
        "chung_buessem_K": float(cb_K)
    }
