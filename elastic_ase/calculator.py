import os
import sys
import numpy as np
import pandas as pd

from ase import Atoms
from elastic_ase.deformations import get_deformation_specs, build_strain_matrix
from elastic_ase.solvers import solve_stiffness, check_mechanical_stability
from elastic_ase.average import calculate_polycrystalline_averages

EV_TO_GPA = 160.21766208

def check_symmetry(atoms: Atoms, symprec: float = 1e-3, verbose: bool = True):
    """
    Analyzes crystal symmetry of the structure using spglib.
    Returns a dict with space group info, Laue class, independent constants, and recommended workflow.
    """
    try:
        import spglib
        cell_data = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
        dataset = spglib.get_symmetry_dataset(cell_data, symprec=symprec)
        if dataset is None:
            return {"space_group": "Unknown", "number": 0, "recommended": "cubic"}

        spg = dataset.number
        if 195 <= spg <= 230:
            system, laue, recommended = "cubic", "m-3m", "cubic"
            constants = ("C11", "C12", "C44")
        elif 168 <= spg <= 194:
            system, laue, recommended = "hexagonal", "6/mmm", "hexagonal"
            constants = ("C11", "C12", "C13", "C33", "C44")
        elif 143 <= spg <= 167:
            system = "trigonal II" if spg <= 148 else "trigonal I"
            laue = "-3" if spg <= 148 else "-3m"
            recommended = "trigonal_2"
            constants = ("C11", "C12", "C13", "C33", "C44", "C14", "C15")
        elif 89 <= spg <= 142:
            system, laue, recommended = "tetragonal I", "4/mmm", "tetragonal_1"
            constants = ("C11", "C12", "C13", "C33", "C44", "C66")
        elif 75 <= spg <= 88:
            system, laue, recommended = "tetragonal II", "4/m", "tetragonal_2"
            constants = ("C11", "C12", "C13", "C33", "C44", "C66", "C16")
        elif 16 <= spg <= 74:
            system, laue, recommended = "orthorhombic", "mmm", "orthogonal"
            constants = ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66")
        elif 3 <= spg <= 15:
            system, laue, recommended = "monoclinic", "2/m", "monoclinic"
            constants = ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66", "C15", "C25", "C35", "C46")
        else:
            system, laue, recommended = "triclinic", "-1", "triclinic"
            constants = tuple(f"C{i}{j}" for i in range(1, 7) for j in range(i, 7))

        report = {
            "formula": atoms.get_chemical_formula(),
            "natoms": len(atoms),
            "space_group": dataset.international,
            "number": spg,
            "crystal_system": system,
            "laue_class": laue,
            "constants": constants,
            "n_constants": len(constants),
            "recommended": recommended,
        }

        if verbose:
            print("=" * 70)
            print(f" ElasticASE Symmetry Check: {report['formula']} ({report['natoms']} atoms)")
            print("=" * 70)
            print(f" Space group       : {report['space_group']} (#{spg}) – {system} (Laue: {laue})")
            print(f" Elastic constants : {len(constants)} independent ({', '.join(constants)})")
            print(f" Recommended mode  : {recommended}")
            print("=" * 70)

        return report
    except Exception as e:
        if verbose:
            print(f"[ElasticASE] Symmetry check failed: {e}")
        return {"space_group": "Unknown", "number": 0, "recommended": "cubic"}

def auto_detect_symmetry(atoms):
    """
    Auto-detect crystal symmetry family using spglib if available.
    """
    try:
        import spglib
        dataset = spglib.get_symmetry_dataset(
            (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
        )
        spg = dataset.number
        if 195 <= spg <= 230:
            return "cubic"
        elif 168 <= spg <= 194:
            return "hexagonal"
        elif 143 <= spg <= 167:
            # Route all trigonal groups to trigonal_2 because it simultaneously fits C14 and C15.
            # This avoids silent errors in '31m'-type groups (P-31m, P-31c, P31m, etc.)
            # and works for any arbitrary in-plane rotation about z.
            return "trigonal_2"
        elif 89 <= spg <= 142:
            return "tetragonal_1"  # 4/mmm (6 independent constants)
        elif 75 <= spg <= 88:
            return "tetragonal_2"  # 4/m (7 independent constants including C16)
        elif 16 <= spg <= 74:
            return "orthogonal"
        elif 3 <= spg <= 15:
            return "monoclinic"
        else:
            return "triclinic"
    except Exception:
        return "cubic"

def validate_orientation(atoms, symmetry, verbose=True):
    """
    Checks for common crystallographic orientation pitfalls:
    1. Trigonal/Hexagonal: c-axis must be along z (cell[2] || z).
    2. Rhombohedral: 3-fold axis along [111] must be in hexagonal setting (c || z).
    3. Monoclinic: unique axis must be b (beta != 90, alpha = gamma = 90).
    """
    warnings = []
    cell = atoms.get_cell()
    angles = cell.angles()
    lengths = cell.lengths()
    sym = symmetry.lower().strip()

    if sym in ["hexagonal", "hex", "trigonal_1", "trigonal_2", "trigonal1", "trigonal2", "tetragonal_1", "tetragonal_2", "tetragonal"]:
        # Rhombohedral check: a=b=c, alpha=beta=gamma != 90
        if sym.startswith("trigonal") and abs(lengths[0] - lengths[1]) < 0.05 and abs(lengths[1] - lengths[2]) < 0.05:
            if abs(angles[0] - angles[1]) < 0.5 and abs(angles[1] - angles[2]) < 0.5 and abs(angles[0] - 90.0) > 1.0:
                warnings.append(
                    f"Rhombohedral cell detected (angles = {angles[0]:.1f}°). "
                    "Trigonal calculations require the 3-fold axis along Cartesian z. "
                    "Please convert to conventional hexagonal setting (c || z)."
                )

        # c-axis along z check
        c_vec = cell[2]
        c_len = np.linalg.norm(c_vec)
        if c_len > 1e-4:
            cos_z = abs(c_vec[2]) / c_len
            if cos_z < 0.90:
                warnings.append(
                    f"{symmetry.capitalize()} c-axis is tilted (|cos theta_z| = {cos_z:.2f} < 0.90). "
                    "The c-axis must be oriented parallel to Cartesian z."
                )

    elif sym in ["monoclinic"]:
        alpha, beta, gamma = angles
        if abs(gamma - 90.0) > 1.5 and abs(beta - 90.0) < 1.0:
            warnings.append(
                f"Monoclinic cell has unique axis c (gamma = {gamma:.1f}°, beta = {beta:.1f}°). "
                "ElasticASE assumes standard unique axis b (alpha=gamma=90°, beta != 90°). "
                "Please permute axes so the 2-fold unique axis is along y."
            )
        elif abs(alpha - 90.0) > 1.5 and abs(beta - 90.0) < 1.0:
            warnings.append(
                f"Monoclinic cell has unique axis a (alpha = {alpha:.1f}°). "
                "ElasticASE assumes standard unique axis b. Permute axes so unique axis is along y."
            )

    if verbose:
        for w in warnings:
            print(f"[ElasticASE Warning] {w}")

    return warnings

def standardize_and_reorient(atoms: Atoms, symprec: float = 1e-3, verbose: bool = True):
    """
    Analyzes cell orientation and automatically standardizes/reorients the structure
    using spglib if a non-standard crystallographic orientation is detected:
    - Rhombohedral cells (3-fold along [111]) -> conventional hexagonal setting (c || z).
    - Tilted hexagonal/trigonal cells (c tilted or along x/y) -> standard setting (c || z, a || x).
    - Monoclinic cells with unique axis c or a -> standard setting (unique axis b || y).
    - Rotated cubic/orthorhombic cells -> standard Cartesian alignment.
    """
    try:
        import spglib
        cell_data = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
        dataset = spglib.get_symmetry_dataset(cell_data, symprec=symprec)
        if dataset is None:
            return atoms, False

        spg = dataset.number
        std_cell_data = spglib.standardize_cell(cell_data, to_primitive=False, symprec=symprec)
        if std_cell_data is None:
            return atoms, False

        std_cell = std_cell_data[0]
        # Check if cell changed substantially in orientation or basis
        if not np.allclose(atoms.get_cell(), std_cell, atol=1e-3):
            new_atoms = Atoms(
                numbers=std_cell_data[2],
                scaled_positions=std_cell_data[1],
                cell=std_cell,
                pbc=True
            )
            new_atoms.calc = atoms.calc
            if verbose:
                print(f"[ElasticASE] Reorientation: Standardized cell to conventional setting (Space group: {dataset.international} #{spg})")
            return new_atoms, True

        return atoms, False
    except Exception as e:
        if verbose:
            print(f"[ElasticASE] Notice: auto-standardization skipped ({e})")
        return atoms, False

class ElasticASE:
    """
    ElasticASE: Compute second-order elastic constants for ASE Atoms objects with ANY ASE calculator.
    Direct ASE successor to elastic_vasp, operating 100% in-memory with analytical accuracy.
    """
    def __init__(
        self,
        atoms: Atoms,
        symmetry: str = "auto",
        strains = None,
        optimizer = None,
        fmax: float = 1e-3,
        reorient: bool = True,
        symprec: float = 1e-3,
        verbose: bool = True
    ):
        if atoms.calc is None:
            raise ValueError("The input ASE Atoms object must have a calculator attached (atoms.calc = ...).")

        self.verbose = verbose
        self.optimizer = optimizer
        self.fmax = fmax
        self.symprec = symprec
        self.reoriented = False

        self.atoms = atoms.copy()
        self.atoms.calc = atoms.calc

        # Automatically decide and reorient structure if requested
        if reorient:
            self.atoms, self.reoriented = standardize_and_reorient(
                self.atoms, symprec=self.symprec, verbose=self.verbose
            )

        if symmetry == "auto":
            self.symmetry = auto_detect_symmetry(self.atoms)
            if self.verbose:
                print(f"[ElasticASE] Auto-detected crystal symmetry: '{self.symmetry}'")
        else:
            self.symmetry = symmetry.lower().strip()

        # Validate crystallographic orientation against common pitfalls
        self.orientation_warnings = validate_orientation(self.atoms, self.symmetry, verbose=self.verbose)

        if strains is None:
            self.strains = np.array([-0.01, -0.005, 0.0, 0.005, 0.01])
        else:
            self.strains = np.array(sorted(strains))

        self.specs = get_deformation_specs(self.symmetry)
        self.V0 = self.atoms.get_volume()
        self.E0 = None
        self.coeffs = {}
        self.r2_scores = {}
        self.stiffness = None
        self.averages = None
        self.directional = None
        self.is_stable = False
        self.stability_criteria = ""

    def calculate(self):
        """
        Executes strain-energy evaluations, polynomial fitting, stiffness assembly, and averaging.
        """
        if self.verbose:
            print("=" * 70)
            print(f" ElasticASE: Elastic Constants Calculation ({self.symmetry.upper()})")
            print("=" * 70)
            print(f"Initial Cell Volume V0 : {self.V0:.4f} Å³ ({self.V0 / len(self.atoms):.4f} Å³/atom)")
            print(f"Strain values          : {self.strains.tolist()}")
            print(f"Deformation modes      : {len(self.specs)}")
            print("-" * 70)

        # Baseline energy at zero strain
        self.E0 = self.atoms.get_potential_energy()

        # Iterate over all deformation modes
        for mode_name, def_func in self.specs:
            energies = []
            for s in self.strains:
                d1, d2, d3, d4, d5, d6 = def_func(float(s))
                strain_mat = build_strain_matrix(d1, d2, d3, d4, d5, d6)

                # Deform cell
                deformed = self.atoms.copy()
                deformed.calc = self.atoms.calc
                orig_cell = self.atoms.get_cell()
                deformed.set_cell(orig_cell @ strain_mat.T, scale_atoms=True)

                # Optional internal coordinate relaxation (for non-affine multi-atom unit cells)
                if self.optimizer is not None and abs(s) > 1e-6:
                    opt = self.optimizer(deformed, logfile=None)
                    opt.run(fmax=self.fmax)

                E = deformed.get_potential_energy()
                energies.append(E)

            energies = np.array(energies)
            delta_E = energies - self.E0

            # Quadratic fit: delta_E / V0 = a * s^2 + b * s + c
            y = delta_E / self.V0
            poly_coeffs = np.polyfit(self.strains, y, 2)
            a = poly_coeffs[0]
            b = poly_coeffs[1]

            # In elastic_vasp convention: coeff = a * 160.21766208
            coeff_val = float(a * EV_TO_GPA)
            self.coeffs[mode_name] = coeff_val

            # Residual stress warning if |b| is non-trivial (> 0.2 GPa)
            residual_stress = abs(float(b * EV_TO_GPA))
            if residual_stress > 0.2 and self.verbose:
                print(f"  [Notice] Mode {mode_name}: residual stress |sigma_0| = {residual_stress:.2f} GPa")

            # Compute R^2 goodness of fit
            poly = np.poly1d(poly_coeffs)
            y_pred = poly(self.strains)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            ss_res = np.sum((y - y_pred) ** 2)
            r2 = float(1.0 - (ss_res / ss_tot)) if ss_tot > 1e-15 else 1.0
            self.r2_scores[mode_name] = r2

            if self.verbose:
                print(f"  Mode {mode_name:<18} -> Coeff = {coeff_val:10.2f} GPa  (R² = {r2:.6f})")

        # Assemble 6x6 stiffness matrix
        self.stiffness = solve_stiffness(self.symmetry, self.coeffs)
        self.is_stable, self.stability_criteria = check_mechanical_stability(self.stiffness, self.symmetry)
        self.averages = calculate_polycrystalline_averages(self.stiffness)
        self.directional = self.directional_analysis()

        if self.verbose:
            self._print_report()

        return {
            "stiffness_matrix": self.stiffness,
            "is_stable": self.is_stable,
            "stability_criteria": self.stability_criteria,
            "polycrystalline_averages": self.averages,
            "directional_properties": self.directional,
            "coefficients": self.coeffs,
            "r2_scores": self.r2_scores
        }

    def _print_report(self):
        print("\n" + "=" * 70)
        print(" 6x6 Stiffness Matrix C_ij (GPa):")
        print("=" * 70)
        df_C = pd.DataFrame(
            np.round(self.stiffness, 2),
            index=[f"C{i}" for i in range(1, 7)],
            columns=[f"C{j}" for j in range(1, 7)]
        )
        print(df_C.to_string())

        print("\n" + "-" * 70)
        print(" Born Mechanical Stability:")
        print(f"  Criteria : {self.stability_criteria}")
        if self.is_stable:
            print("  Status   : \033[1;32mMECHANICALLY STABLE\033[0m")
        else:
            print("  Status   : \033[1;31mMECHANICALLY UNSTABLE\033[0m")

        print("\n" + "-" * 70)
        print(" Polycrystalline Moduli (Voigt-Reuss-Hill):")
        print(f"  Bulk Modulus K_VRH   : {self.averages['K_VRH']:8.2f} GPa  (Voigt: {self.averages['K_voigt']:.2f}, Reuss: {self.averages['K_reuss']:.2f})")
        print(f"  Shear Modulus G_VRH  : {self.averages['G_VRH']:8.2f} GPa  (Voigt: {self.averages['G_voigt']:.2f}, Reuss: {self.averages['G_reuss']:.2f})")
        print(f"  Young's Modulus E    : {self.averages['E_VRH']:8.2f} GPa")
        print(f"  Poisson's Ratio nu   : {self.averages['nu_VRH']:8.3f}")
        print(f"  Pugh Ratio (K/G)     : {self.averages['Pugh_ratio']:8.2f}  (> 1.75 indicates ductile)")
        print(f"  Cauchy Pressure      : {self.averages['Cauchy_pressure']:8.2f} GPa")
        print(f"  Zener Anisotropy A   : {self.averages['Zener_anisotropy']:8.2f}")
        print(f"  Universal Anisotropy : A^U = {self.averages['universal_anisotropy']:6.3f}  (0 indicates isotropic)")

        if self.directional:
            print("\n" + "-" * 70)
            print(" Directional Elastic Extremes & Anisotropy (elastic-vis):")
            print(f"  Young's Modulus E    : min = {self.directional['E_min']:8.2f} GPa, max = {self.directional['E_max']:8.2f} GPa  (Anisotropy: {self.directional['E_anisotropy']:.2f})")
            print(f"  Shear Modulus G      : min = {self.directional['G_min']:8.2f} GPa, max = {self.directional['G_max']:8.2f} GPa  (Anisotropy: {self.directional['G_anisotropy']:.2f})")
            aux = "YES (Auxetic behavior)" if self.directional["is_auxetic"] else "No"
            print(f"  Poisson's Ratio nu   : min = {self.directional['nu_min']:8.3f},     max = {self.directional['nu_max']:8.3f}      (Negative nu: {aux})")
            neg_lc = "YES (Negative LC)" if self.directional["is_negative_LC"] else "No"
            print(f"  Linear Compress. LC  : min = {self.directional['LC_min']:8.2f},     max = {self.directional['LC_max']:8.2f} TPa⁻¹ (Negative LC: {neg_lc})")
        print("=" * 70)

    @property
    def compliance(self):
        """6x6 Voigt compliance matrix S = C^-1 in 1/GPa."""
        if self.stiffness is None:
            raise RuntimeError("Must call calculate() before accessing compliance matrix.")
        return np.linalg.inv(self.stiffness)

    def directional_analysis(self):
        """
        Calculates directional elastic properties and extremes without requiring any plotting:
        - Young's modulus E: min, max, anisotropy (E_max / E_min)
        - Shear modulus G: min, max, anisotropy (G_max / G_min)
        - Poisson's ratio nu: min, max, auxetic check
        - Linear compressibility LC: min, max, negative LC check
        """
        if self.stiffness is None:
            raise RuntimeError("Must call calculate() before performing directional analysis.")

        try:
            material = self.to_elastic_vis()
            from elastic_vis.core import minimize, maximize

            _, Emin = minimize(lambda x: material.Young(x), 2)
            _, Emax = maximize(lambda x: material.Young(x), 2)
            E_aniso = float(Emax / Emin) if Emin > 1e-6 else float("nan")

            _, Gmin = minimize(lambda x: material.shear(x), 3)
            _, Gmax = maximize(lambda x: material.shear(x), 3)
            G_aniso = float(Gmax / Gmin) if Gmin > 1e-6 else float("nan")

            _, nu_min = minimize(lambda x: material.Poisson(x), 3)
            _, nu_max = maximize(lambda x: material.Poisson(x), 3)

            _, lc_min = minimize(lambda x: material.LC(x), 2)
            _, lc_max = maximize(lambda x: material.LC(x), 2)

            self.directional = {
                "E_min": float(Emin),
                "E_max": float(Emax),
                "E_anisotropy": float(E_aniso),
                "G_min": float(Gmin),
                "G_max": float(Gmax),
                "G_anisotropy": float(G_aniso),
                "nu_min": float(nu_min),
                "nu_max": float(nu_max),
                "is_auxetic": bool(nu_min < 0.0),
                "LC_min": float(lc_min),
                "LC_max": float(lc_max),
                "is_negative_LC": bool(lc_min < 0.0),
            }
        except Exception:
            self.directional = None

        return self.directional

    def Young(self, theta: float, phi: float) -> float:
        """Directional Young's modulus along spherical angles (theta, phi) in GPa."""
        return self.to_elastic_vis().Young([theta, phi])

    def linear_compressibility(self, theta: float, phi: float) -> float:
        """Directional Linear Compressibility along spherical angles (theta, phi) in 1000/GPa (TPa^-1)."""
        return self.to_elastic_vis().LC([theta, phi])

    def shear(self, theta: float, phi: float) -> tuple:
        """Extremes of Shear Modulus (min, max) along direction (theta, phi) in GPa."""
        return self.to_elastic_vis().shear2D([theta, phi])

    def poisson(self, theta: float, phi: float) -> tuple:
        """Extremes of Poisson's Ratio (min, max) along direction (theta, phi)."""
        return self.to_elastic_vis().Poisson2D([theta, phi])

    def to_elastic_vis(self):
        """
        Converts the calculated 6x6 stiffness matrix to an elastic_vis.Elastic object.
        """
        if self.stiffness is None:
            raise RuntimeError("Must call calculate() before converting to elastic_vis.")
        from elastic_vis import Elastic
        return Elastic(self.stiffness)

    def visualize(
        self,
        properties=("Young", "Shear"),
        dimensions=("2D", "3D"),
        prefix="elastic",
        style_config=None,
        show=False,
        output_format="png"
    ):
        """
        Visualizes directional elastic properties using elastic-vis:
        - properties: list/tuple containing 'Young', 'Shear', 'Poisson', 'LC' (or 'all')
        - dimensions: list/tuple containing '2D', '3D' (or 'both')
        - prefix: filename prefix for saved plots
        - style_config: optional dict or path to style JSON file
        - show: whether to display plots interactively
        - output_format: 'png', 'pdf', or 'svg'
        """
        if self.stiffness is None:
            raise RuntimeError("Must call calculate() before visualising results.")

        import json
        from elastic_vis import (
            Elastic,
            plot_youngs_2d,
            plot_youngs_3d,
            plot_shear_2d,
            plot_shear_3d,
            plot_poisson_2d,
            plot_poisson_3d,
            plot_lc_2d,
            plot_lc_3d
        )

        if isinstance(style_config, str):
            with open(style_config, "r") as f:
                style_config = json.load(f)

        material = Elastic(self.stiffness)

        if isinstance(properties, str):
            if properties.lower() == "all":
                properties = ["Young", "Shear", "Poisson", "LC"]
            else:
                properties = [properties]
        prop_set = {p.lower() for p in properties}

        if isinstance(dimensions, str):
            if dimensions.lower() in ["both", "all"]:
                dimensions = ["2D", "3D"]
            else:
                dimensions = [dimensions]
        dim_set = {d.upper() for d in dimensions}

        generated_files = []

        if "young" in prop_set:
            if "2d" in [d.lower() for d in dim_set]:
                plot_youngs_2d(material, prefix=prefix, style_config=style_config, show=show, output_format=output_format)
                generated_files.append(f"{prefix}_Youngs_Modulus_Polar_2D.{output_format}")
            if "3d" in [d.lower() for d in dim_set]:
                plot_youngs_3d(material, prefix=prefix, style_config=style_config, show=show, output_format=output_format)
                generated_files.append(f"{prefix}_Youngs_Modulus_3D.{output_format}")

        if "shear" in prop_set:
            if "2d" in [d.lower() for d in dim_set]:
                plot_shear_2d(material, prefix=prefix, style_config=style_config, show=show, output_format=output_format)
                generated_files.append(f"{prefix}_Shear_Modulus_Polar_2D.{output_format}")
            if "3d" in [d.lower() for d in dim_set]:
                plot_shear_3d(material, prefix=prefix, style_config=style_config, show=show, output_format=output_format)
                generated_files.append(f"{prefix}_Shear_Modulus_3D.{output_format}")

        if "poisson" in prop_set:
            if "2d" in [d.lower() for d in dim_set]:
                plot_poisson_2d(material, prefix=prefix, style_config=style_config, show=show, output_format=output_format)
                generated_files.append(f"{prefix}_Poisson_Ratio_Polar_2D.{output_format}")
            if "3d" in [d.lower() for d in dim_set]:
                plot_poisson_3d(material, prefix=prefix, style_config=style_config, show=show, output_format=output_format)
                generated_files.append(f"{prefix}_Poisson_Ratio_3D.{output_format}")

        if "lc" in prop_set or "linear_compressibility" in prop_set:
            if "2d" in [d.lower() for d in dim_set]:
                plot_lc_2d(material, prefix=prefix, style_config=style_config, show=show, output_format=output_format)
                generated_files.append(f"{prefix}_Linear_Compressibility_Polar_2D.{output_format}")
            if "3d" in [d.lower() for d in dim_set]:
                plot_lc_3d(material, prefix=prefix, style_config=style_config, show=show, output_format=output_format)
                generated_files.append(f"{prefix}_Linear_Compressibility_3D.{output_format}")

        if self.verbose:
            print(f"[ElasticASE] Visualization saved: {', '.join(generated_files)}")

        return generated_files

    def plot_youngs_3d(self, prefix="elastic", **kwargs):
        return self.visualize(properties=["Young"], dimensions=["3D"], prefix=prefix, **kwargs)

    def plot_youngs_2d(self, prefix="elastic", **kwargs):
        return self.visualize(properties=["Young"], dimensions=["2D"], prefix=prefix, **kwargs)

    def plot_shear_3d(self, prefix="elastic", **kwargs):
        return self.visualize(properties=["Shear"], dimensions=["3D"], prefix=prefix, **kwargs)

    def plot_shear_2d(self, prefix="elastic", **kwargs):
        return self.visualize(properties=["Shear"], dimensions=["2D"], prefix=prefix, **kwargs)

    def plot_poisson_3d(self, prefix="elastic", **kwargs):
        return self.visualize(properties=["Poisson"], dimensions=["3D"], prefix=prefix, **kwargs)

    def plot_poisson_2d(self, prefix="elastic", **kwargs):
        return self.visualize(properties=["Poisson"], dimensions=["2D"], prefix=prefix, **kwargs)

    def plot_lc_3d(self, prefix="elastic", **kwargs):
        return self.visualize(properties=["LC"], dimensions=["3D"], prefix=prefix, **kwargs)

    def plot_lc_2d(self, prefix="elastic", **kwargs):
        return self.visualize(properties=["LC"], dimensions=["2D"], prefix=prefix, **kwargs)

    def export_matrix(self, filepath: str):
        """
        Exports the 6x6 stiffness matrix as a plain space-separated text file,
        compatible with elastic-vis and numpy.loadtxt.
        """
        if self.stiffness is None:
            raise RuntimeError("Must call calculate() before exporting matrix.")
        np.savetxt(filepath, self.stiffness, fmt="%.6f")
        if self.verbose:
            print(f"[ElasticASE] Exported 6x6 stiffness matrix to {filepath}")

    def __repr__(self):
        status = "calculated" if self.stiffness is not None else "uncalculated"
        return f"<ElasticASE(symmetry='{self.symmetry}', status='{status}')>"

    def __str__(self):
        if self.stiffness is None:
            return f"ElasticASE(symmetry='{self.symmetry}', status='uncalculated')"
        c11 = self.stiffness[0, 0]
        stable = "Stable" if self.is_stable else "Unstable"
        return f"ElasticASE(symmetry='{self.symmetry}', C11={c11:.2f} GPa, {stable})"
