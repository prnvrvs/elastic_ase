"""From energies to elastic constants: pure functions (no files, no printing).

For every deformation the energy per volume is fitted with (E - E0)/V0 = a d^2 + b d + c; the
curvatures a are solved for the independent constants of the table (see the theory page of
the documentation).
"""
import numpy as np

from elastic_ase import auto
from elastic_ase.average import calculate_polycrystalline_averages
from elastic_ase.properties import derived_properties
from elastic_ase.results import ElasticResults
from elastic_ase.rundir import GROUND_STATE, strained_key
from elastic_ase.solvers import check_mechanical_stability

EV_TO_GPA = 160.21766208
RESIDUAL_STRESS_TOL = 0.5  # GPa: warn above this residual stress of the reference cell


def fit_curves(system, strains, energies, volume):
    """Quadratic fit of every energy-strain curve.

    Returns (coefficients, r2_scores, fits, linear_terms): the curvature a (GPa) and R^2 per
    deformation, the data and fitted polynomial per deformation (for plots), and the
    (strain direction e, linear term b in GPa) pairs used for the residual stress.
    """
    strains = np.asarray(strains, float)
    e0 = energies[GROUND_STATE]
    coefficients, r2_scores, fits, linear_terms = {}, {}, {}, []
    for deformation in system.deformations:
        name = deformation.folder
        energy = np.array([e0 if s == 0.0 else energies[strained_key(name, float(s))] for s in strains])
        y = (energy - e0) / volume
        poly = np.polyfit(strains, y, 2)
        coefficients[name] = float(poly[0] * EV_TO_GPA)

        residual = y - np.polyval(poly, strains)
        ss_tot = np.sum((y - y.mean()) ** 2)
        r2_scores[name] = float(1.0 - np.sum(residual ** 2) / ss_tot) if ss_tot > 1e-15 else 1.0
        fits[name] = {"strains": strains.tolist(), "y_GPa": (y * EV_TO_GPA).tolist(),
                      "coefficients_GPa": (poly * EV_TO_GPA).tolist(), "r2": r2_scores[name]}

        # Linear term = sigma_0 . e. From a cubic fit when possible: in a quadratic fit the
        # anharmonic cubic term leaks into it (about 0.3 GPa of fake pressure for Al at +-3 %).
        b = np.polyfit(strains, y, 3)[2] if len(np.unique(strains)) >= 4 else poly[1]
        linear_terms.append((np.array(deformation.strain, float), float(b * EV_TO_GPA)))
    return coefficients, r2_scores, fits, linear_terms


def residual_stress(system, linear_terms):
    """Residual stress of the reference cell from the linear terms of the fits.

    b_j = sigma_0 . e_j for every deformation j; solving gives the parts of sigma_0 the
    deformations determine. The pressure P = -tr(sigma_0)/3 is exact whenever (1,1,1,0,0,0) lies
    in the span of the e_j (true for every named table), else None.
    """
    e = np.array([t[0] for t in linear_terms])
    b = np.array([t[1] for t in linear_terms])
    sigma, *_ = np.linalg.lstsq(e, b, rcond=None)
    hydro = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    coef, *_ = np.linalg.lstsq(e.T, hydro, rcond=None)
    pressure = float(-hydro @ sigma / 3.0) if np.abs(e.T @ coef - hydro).max() < 1e-8 else None
    along = np.abs(b) / np.linalg.norm(e, axis=1)
    return {
        "pressure_GPa": pressure,  # positive = compressed reference cell
        "max_along_deformation_GPa": float(along.max()),
        "per_deformation_GPa": {d.folder: float(v) for d, v in zip(system.deformations, along)},
    }


def residual_stress_message(stress, tolerance=RESIDUAL_STRESS_TOL):
    """Warning text if the reference cell is noticeably stressed, else None."""
    worst = stress["max_along_deformation_GPa"]
    if worst <= tolerance:
        return None
    pressure = stress["pressure_GPa"]
    text = f"pressure {pressure:+.2f} GPa" if pressure is not None else f"up to {worst:.2f} GPa"
    return (f"The reference cell is not stress-free ({text}, from the linear terms of the "
            "energy-strain fits). The elastic constants are biased by a similar amount: relax "
            "the structure (cell and atoms) with the same calculator and settings first.")


def standard_frame(stiffness, rotation, space_group_number, constants):
    """(C in the standard setting, independent constants of the space group there)."""
    if rotation is None:
        return stiffness.copy(), dict(constants)
    c_standard = auto.to_standard(stiffness, rotation)
    standard_constants = dict(constants)
    if space_group_number > 0:
        try:
            _, standard_constants, _ = auto.standard_constants(c_standard, space_group_number)
        except Exception:
            standard_constants = dict(constants)
    return c_standard, standard_constants


def analyse(system, strains, energies, atoms, symmetry, standard_rotation=None,
            space_group_number=0, frame="cartesian"):
    """ElasticResults from the energies of all runs (``energies[run name]`` in eV).

    ``atoms`` is the reference structure (its volume and masses are used), ``symmetry`` the
    table name (for the stability description and the Kleinman parameter). Directional
    properties (elastic-vis) are not computed here.
    """
    coefficients, r2_scores, fits, linear_terms = fit_curves(system, strains, energies, atoms.get_volume())
    constants, stiffness = system.solve(coefficients)
    c_standard, standard_constants = standard_frame(stiffness, standard_rotation, space_group_number, constants)

    # Stability and the Zener/Cauchy values read entries such as C12 and C44, which have their
    # usual meaning only with the crystal axes along x, y, z: use the standard setting.
    is_stable, criteria = check_mechanical_stability(c_standard, symmetry)
    averages = calculate_polycrystalline_averages(c_standard)
    return ElasticResults(
        stiffness_matrix=c_standard if frame == "standard" else stiffness,
        stiffness_cartesian=stiffness,
        stiffness_standard=c_standard,
        standard_constants=standard_constants,
        standard_rotation=standard_rotation,
        independent_constants=constants,
        is_stable=is_stable,
        stability_criteria=criteria,
        polycrystalline_averages=averages,
        directional_properties=None,
        derived_properties=derived_properties(averages, atoms, c_standard, symmetry),
        residual_stress=residual_stress(system, linear_terms),
        coefficients=coefficients,
        r2_scores=r2_scores,
        energies=dict(energies),
        fits=fits,
    )
