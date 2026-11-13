"""Properties derived from the elastic constants and the structure.

All moduli in GPa, velocities in m/s, temperatures in K.

* Sound velocities from the Hill moduli (isotropic polycrystal):
  v_l = sqrt((K + 4G/3) / rho), v_t = sqrt(G / rho), mean v_m = [(2/v_t^3 + 1/v_l^3) / 3]^(-1/3).
* Debye temperature (Anderson, J. Phys. Chem. Solids 24, 909 (1963)):
  theta_D = h/k_B * (3 n / (4 pi V))^(1/3) * v_m, with n/V the number of atoms per volume.
* Kleinman parameter (cubic only): zeta = (C11 + 8 C12) / (7 C11 + 2 C12).
* Empirical Vickers hardness, with k = G/K (Pugh ratio inverse):
  Chen et al., Intermetallics 19, 1275 (2011): H_v = 2 (k^2 G)^0.585 - 3;
  Tian et al., Int. J. Refract. Met. Hard Mater. 33, 93 (2012): H_v = 0.92 k^1.137 G^0.708.
  Both were fitted to hard covalent/ionic solids and are rough estimates; for soft metals
  the Chen value can even be negative.
"""
import numpy as np

AMU_PER_A3_TO_KG_PER_M3 = 1.66053906660e3
H_OVER_KB = 6.62607015e-34 / 1.380649e-23  # K s


def derived_properties(averages, atoms, stiffness_standard=None, symmetry=None):
    """Dictionary of density, sound velocities, Debye temperature, hardness estimates, and
    for cubic crystals the Kleinman parameter."""
    volume = atoms.get_volume()  # A^3
    rho = atoms.get_masses().sum() / volume * AMU_PER_A3_TO_KG_PER_M3  # kg/m^3
    K, G = averages["K_VRH"] * 1e9, averages["G_VRH"] * 1e9  # Pa

    out = {"density_kg_m3": float(rho)}
    if G > 0 and K + 4.0 * G / 3.0 > 0:
        v_l = np.sqrt((K + 4.0 * G / 3.0) / rho)
        v_t = np.sqrt(G / rho)
        v_m = ((2.0 / v_t ** 3 + 1.0 / v_l ** 3) / 3.0) ** (-1.0 / 3.0)
        n_per_m3 = len(atoms) / (volume * 1e-30)
        theta = H_OVER_KB * (3.0 * n_per_m3 / (4.0 * np.pi)) ** (1.0 / 3.0) * v_m
        out.update(v_longitudinal_m_s=float(v_l), v_transverse_m_s=float(v_t),
                   v_mean_m_s=float(v_m), debye_temperature_K=float(theta))
    else:
        out.update(v_longitudinal_m_s=float("nan"), v_transverse_m_s=float("nan"),
                   v_mean_m_s=float("nan"), debye_temperature_K=float("nan"))

    Kg, Gg = averages["K_VRH"], averages["G_VRH"]
    if Kg > 0 and Gg > 0:
        k = Gg / Kg
        out["hardness_chen_GPa"] = float(2.0 * (k * k * Gg) ** 0.585 - 3.0)
        out["hardness_tian_GPa"] = float(0.92 * k ** 1.137 * Gg ** 0.708)
    else:
        out["hardness_chen_GPa"] = out["hardness_tian_GPa"] = float("nan")

    if symmetry == "cubic" and stiffness_standard is not None:
        c11, c12 = stiffness_standard[0, 0], stiffness_standard[0, 1]
        out["kleinman_parameter"] = float((c11 + 8.0 * c12) / (7.0 * c11 + 2.0 * c12))
    return out
