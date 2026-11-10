import os
import sys
import numpy as np
import pandas as pd

from ase import Atoms
from elastic_ase.deformations import get_deformation_specs, build_strain_matrix
from elastic_ase.solvers import solve_stiffness, check_mechanical_stability
from elastic_ase.average import calculate_polycrystalline_averages

EV_TO_GPA = 160.21766208

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
            return "trigonal_1"
        elif 75 <= spg <= 142:
            return "orthogonal"  # tetragonal treated as orthogonal
        elif 16 <= spg <= 74:
            return "orthogonal"
        elif 3 <= spg <= 15:
            return "monoclinic"
        else:
            return "triclinic"
    except Exception:
        return "cubic"

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
        verbose: bool = True
    ):
        if atoms.calc is None:
            raise ValueError("The input ASE Atoms object must have a calculator attached (atoms.calc = ...).")

        self.atoms = atoms.copy()
        self.atoms.calc = atoms.calc
        self.verbose = verbose
        self.optimizer = optimizer
        self.fmax = fmax

        if symmetry == "auto":
            self.symmetry = auto_detect_symmetry(atoms)
            if self.verbose:
                print(f"[ElasticASE] Auto-detected crystal symmetry: '{self.symmetry}'")
        else:
            self.symmetry = symmetry.lower().strip()

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

        if self.verbose:
            self._print_report()

        return {
            "stiffness_matrix": self.stiffness,
            "is_stable": self.is_stable,
            "stability_criteria": self.stability_criteria,
            "polycrystalline_averages": self.averages,
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
        print("=" * 70)

    def __repr__(self):
        status = "calculated" if self.stiffness is not None else "uncalculated"
        return f"<ElasticASE(symmetry='{self.symmetry}', status='{status}')>"

    def __str__(self):
        if self.stiffness is None:
            return f"ElasticASE(symmetry='{self.symmetry}', status='uncalculated')"
        c11 = self.stiffness[0, 0]
        stable = "Stable" if self.is_stable else "Unstable"
        return f"ElasticASE(symmetry='{self.symmetry}', C11={c11:.2f} GPa, {stable})"
