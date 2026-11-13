"""Plain-text report of an analysis."""
import numpy as np


def format_matrix(matrix, decimals=2):
    """6x6 matrix as text with C1..C6 row and column labels."""
    labels = [f"C{i}" for i in range(1, 7)]
    cells = [[f"{v:.{decimals}f}" for v in row] for row in np.asarray(matrix, float)]
    width = max(len(c) for row in cells for c in row)
    lines = ["    " + "  ".join(label.rjust(width) for label in labels)]
    for label, row in zip(labels, cells):
        lines.append(f"{label:<4}" + "  ".join(c.rjust(width) for c in row))
    return "\n".join(lines)


def format_report(results, frame="cartesian", space_group="", color=True):
    """The report printed after an analysis (``results`` is an ElasticResults)."""
    bar, rule = "=" * 70, "-" * 70
    green, red, reset = ("\033[1;32m", "\033[1;31m", "\033[0m") if color else ("", "", "")
    r = results
    out = ["", bar, f" 6x6 Stiffness Matrix C_ij (GPa) [{'Standard Setting' if frame == 'standard' else 'Cartesian Frame'}]:",
           bar, format_matrix(r.stiffness_matrix)]

    if frame != "standard" and not np.allclose(r.stiffness_cartesian, r.stiffness_standard, atol=1e-2):
        out += ["", rule, " Standard Crystallographic Setting (C_std):"]
        if space_group:
            out.append(f"  Space group : {space_group}")
        if r.standard_constants:
            out.append("  Constants   : " + ", ".join(f"{k} = {v:.2f}" for k, v in r.standard_constants.items()))
        out.append(format_matrix(r.stiffness_standard))

    status = f"{green}MECHANICALLY STABLE{reset}" if r.is_stable else f"{red}MECHANICALLY UNSTABLE{reset}"
    out += ["", rule, " Born Mechanical Stability:", f"  Criteria : {r.stability_criteria}", f"  Status   : {status}"]

    a = r.polycrystalline_averages
    out += ["", rule, " Polycrystalline Moduli (Voigt-Reuss-Hill):",
            f"  Bulk Modulus K_VRH   : {a['K_VRH']:8.2f} GPa  (Voigt: {a['K_voigt']:.2f}, Reuss: {a['K_reuss']:.2f})",
            f"  Shear Modulus G_VRH  : {a['G_VRH']:8.2f} GPa  (Voigt: {a['G_voigt']:.2f}, Reuss: {a['G_reuss']:.2f})",
            f"  Young's Modulus E    : {a['E_VRH']:8.2f} GPa",
            f"  Poisson's Ratio nu   : {a['nu_VRH']:8.3f}",
            f"  Pugh Ratio (K/G)     : {a['Pugh_ratio']:8.2f}  (> 1.75 usually ductile)",
            f"  Cauchy Pressure      : {a['Cauchy_pressure']:8.2f} GPa",
            f"  Zener Anisotropy A   : {a['Zener_anisotropy']:8.2f}",
            f"  Universal Anisotropy : A^U = {a['universal_anisotropy']:6.3f}  (0 for isotropic)"]

    p = r.derived_properties
    out += ["", rule, " Derived Properties:",
            f"  Density              : {p['density_kg_m3']:8.1f} kg/m³",
            f"  Sound velocity       : longitudinal {p['v_longitudinal_m_s']:.0f} m/s, transverse "
            f"{p['v_transverse_m_s']:.0f} m/s, mean {p['v_mean_m_s']:.0f} m/s",
            f"  Debye temperature    : {p['debye_temperature_K']:8.1f} K  (Anderson)"]
    if "kleinman_parameter" in p:
        out.append(f"  Kleinman parameter   : {p['kleinman_parameter']:8.3f}")
    out.append(f"  Vickers hardness     : {p['hardness_chen_GPa']:.2f} GPa (Chen), "
               f"{p['hardness_tian_GPa']:.2f} GPa (Tian)  [empirical, rough]")
    if r.residual_stress["pressure_GPa"] is not None:
        out.append(f"  Residual pressure    : {r.residual_stress['pressure_GPa']:+8.2f} GPa  (from the fits; should be ~0)")

    d = r.directional_properties
    if d:
        out += ["", rule, " Directional Elastic Extremes & Anisotropy (elastic-vis):",
                f"  Young's Modulus E    : min = {d['E_min']:8.2f} GPa, max = {d['E_max']:8.2f} GPa  (Anisotropy: {d['E_anisotropy']:.2f})",
                f"  Shear Modulus G      : min = {d['G_min']:8.2f} GPa, max = {d['G_max']:8.2f} GPa  (Anisotropy: {d['G_anisotropy']:.2f})",
                f"  Poisson's Ratio nu   : min = {d['nu_min']:8.3f},     max = {d['nu_max']:8.3f}      "
                f"(Negative nu: {'YES (auxetic)' if d['is_auxetic'] else 'No'})",
                f"  Linear Compress. LC  : min = {d['LC_min']:8.2f},     max = {d['LC_max']:8.2f} TPa⁻¹ "
                f"(Negative LC: {'YES' if d['is_negative_LC'] else 'No'})"]
    out.append(bar)
    return "\n".join(out)
