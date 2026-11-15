"""Residual stress check, derived properties, fit plots, save/load, QE helpers."""
import os
import tempfile
import unittest
import warnings

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.units import GPa

from elastic_ase import ElasticASE, ElasticConstants, calculate_polycrystalline_averages, derived_properties
from elastic_ase.calculator import _qe_converged

A0 = 3.9943  # EMT equilibrium lattice constant of Al


def run(a, **kw):
    atoms = bulk("Al", "fcc", a=a, cubic=True)
    atoms.calc = EMT()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ec = ElasticASE(atoms, verbose=False, **kw)
        ec.calculate()
    stress_warnings = [w for w in caught if "stress-free" in str(w.message)]
    return ec, atoms, stress_warnings


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


class TestFeatures(unittest.TestCase):
    def test_alias(self):
        self.assertIs(ElasticConstants, ElasticASE)

    def test_residual_pressure_from_fits(self):
        ec, atoms, warned = run(4.05)  # EMT: about 1.56 GPa of tension
        true_p = -atoms.get_stress()[:3].mean() / GPa
        self.assertAlmostEqual(ec.residual_stress["pressure_GPa"], true_p, delta=0.01)
        self.assertEqual(len(warned), 1)

        ec, _, warned = run(A0)
        self.assertLess(abs(ec.residual_stress["pressure_GPa"]), 0.01)
        self.assertEqual(warned, [])

    def test_no_fake_pressure_from_anharmonicity(self):
        # large strains on a relaxed cell: the cubic term must not show up as a pressure
        ec, _, warned = run(A0, strains=[-0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03])
        self.assertLess(abs(ec.residual_stress["pressure_GPa"]), 0.05)
        self.assertEqual(warned, [])

    def test_derived_properties_aluminium(self):
        # VASP (PBE) constants of fcc Al from the tutorial; experiment: v_l ~ 6420 m/s, theta_D ~ 428 K
        C = np.zeros((6, 6))
        C[:3, :3] = 63.4
        np.fill_diagonal(C[:3, :3], 106.6)
        C[3, 3] = C[4, 4] = C[5, 5] = 28.9
        p = derived_properties(calculate_polycrystalline_averages(C), bulk("Al", "fcc", a=4.04, cubic=True), C, "cubic")
        self.assertAlmostEqual(p["density_kg_m3"], 2717.9, delta=1.0)
        self.assertAlmostEqual(p["v_longitudinal_m_s"], 6422, delta=5)
        self.assertAlmostEqual(p["v_transverse_m_s"], 3076, delta=5)
        self.assertAlmostEqual(p["debye_temperature_K"], 404.7, delta=1.0)
        self.assertNotIn("kleinman_parameter", p)   # fcc: no internal strain
        p_si = derived_properties(calculate_polycrystalline_averages(C), bulk("Si", "diamond", a=5.43), C, "cubic")
        self.assertAlmostEqual(p_si["kleinman_parameter"], 0.703, delta=0.002)
        p_hex = derived_properties(calculate_polycrystalline_averages(C), bulk("Al", "fcc", a=4.04), C, "hexagonal")
        self.assertNotIn("kleinman_parameter", p_hex)

    def test_results_contain_new_entries(self):
        ec, _, _ = run(A0)
        for key in ("derived_properties", "residual_stress", "energies"):
            self.assertIn(key, ec.results)
        self.assertEqual(len(ec.results["energies"]), len(ec.runs()))

    def test_plot_fits(self):
        ec, _, _ = run(A0)
        with tempfile.TemporaryDirectory() as tmp:
            path = ec.plot_fits(os.path.join(tmp, "fits.png"))
            self.assertGreater(os.path.getsize(path), 1000)
        fresh = ElasticASE(bulk("Al", "fcc", a=A0, cubic=True), verbose=False)
        with self.assertRaises(RuntimeError):
            fresh.plot_fits()

    def test_save_and_load_reproduce_results(self):
        hcp = bulk("Cu", "hcp", a=2.54, c=4.15)
        hcp.rotate(25, "z", rotate_cell=True)
        for kw in ({}, {"reorient": False}):  # named table after reorientation; orientation-aware auto
            atoms = hcp.copy()
            atoms.calc = EMT()
            ec = ElasticASE(atoms, verbose=False, **kw)
            res = ec.calculate()
            with tempfile.TemporaryDirectory() as tmp:
                loaded = ElasticASE.load(ec.save(os.path.join(tmp, "r.json")))
            np.testing.assert_allclose(loaded.results["stiffness_matrix"], res["stiffness_matrix"], atol=1e-10)
            self.assertEqual(loaded.symmetry, ec.symmetry)
            self.assertAlmostEqual(loaded.properties["debye_temperature_K"], ec.properties["debye_temperature_K"])
        with self.assertRaises(RuntimeError):
            ElasticASE(bulk("Al", "fcc", a=A0, cubic=True), verbose=False).save("x.json")

    def test_qe_convergence_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            ok, bad, other = (os.path.join(tmp, n) for n in ("ok.pwo", "bad.pwo", "OUTCAR"))
            _write(ok, "     Program PWSCF v.7.4 starts\n     convergence has been achieved in  8 iterations\n")
            _write(bad, "     Program PWSCF v.7.4 starts\n     convergence NOT achieved after 100 iterations: stopping\n")
            _write(other, "vasp output\n")
            self.assertTrue(_qe_converged(ok))
            self.assertFalse(_qe_converged(bad))
            self.assertIsNone(_qe_converged(other))

    def test_qe_input(self):
        from ase.calculators.espresso import Espresso, EspressoProfile
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):  # no Espresso calculator: pseudopotentials unknown
                ElasticASE(bulk("Al", "fcc", a=A0), verbose=False).qe_input(tmp)
            atoms = bulk("Al", "fcc", a=A0)
            atoms.calc = Espresso(profile=EspressoProfile(command="pw.x", pseudo_dir="pseudo"),
                                  pseudopotentials={"Al": "Al.UPF"},
                                  input_data={"system": {"ecutwfc": 60}}, kpts=(8, 8, 8))
            ec = ElasticASE(atoms, verbose=False)
            folders = ec.qe_input(tmp)
            self.assertEqual(len(folders), len(ec.runs()))
            with open(os.path.join(folders[1], "espresso.pwi")) as f:
                text = f.read()
            self.assertIn("ecutwfc", text)
            self.assertIn("8 8 8", text)
            self.assertEqual(len(ec.qe_status(tmp)["pending"]), len(ec.runs()))


    def test_results_object(self):
        from elastic_ase.results import ElasticResults
        ec, _, _ = run(A0)
        r = ec.results
        self.assertIsInstance(r, ElasticResults)
        self.assertIs(r["stiffness_matrix"], r.stiffness_matrix)  # works as a dictionary too
        self.assertIn("derived_properties", r)
        self.assertEqual(set(dict(r)), set(r.to_dict()))
        with self.assertRaises(KeyError):
            r["no_such_key"]
        # the old attributes are read-only views of the one results object
        self.assertIs(ec.stiffness, r.stiffness_cartesian)
        self.assertIs(ec.averages, r.polycrystalline_averages)
        self.assertEqual(ec.E0, r.energies["ground_state"])
        with self.assertRaises(AttributeError):
            ec.stiffness = None
        fresh = ElasticASE(bulk("Al", "fcc", a=A0, cubic=True), verbose=False)
        self.assertIsNone(fresh.stiffness)
        self.assertIsNone(fresh.results)

    def test_code_adapters(self):
        from elastic_ase import codes
        self.assertIs(codes.get_code("Espresso"), codes.CODES["qe"])
        self.assertIsNone(codes.get_code(None))
        with self.assertRaises(ValueError):
            codes.get_code("castepp")
        with tempfile.TemporaryDirectory() as tmp:
            bad_qe, bad_vasp = os.path.join(tmp, "x.pwo"), os.path.join(tmp, "OUTCAR")
            _write(bad_qe, "Program PWSCF\nconvergence NOT achieved after 100 iterations\n")
            _write(bad_vasp, "aborting loop EDIFF was not reached (unconverged)\n")
            self.assertIn("QE", codes.not_converged(bad_qe))
            self.assertIn("VASP", codes.not_converged(bad_vasp))
            self.assertIsNone(codes.not_converged(bad_qe, codes.CODES["vasp"]))  # only the given code's check
            # an unconverged pw.x output that ASE cannot parse is still named as unconverged
            from elastic_ase.rundir import read_outputs
            os.makedirs(os.path.join(tmp, "ground_state"))
            _write(os.path.join(tmp, "ground_state", "espresso.pwo"),
                   "Program PWSCF v.7.6 starts\nconvergence NOT achieved after 1 iterations: stopping\n")
            energies, problems = read_outputs(tmp, [("ground_state", bulk("Al"))], "espresso.pwo",
                                              codes.CODES["qe"])
            self.assertEqual(energies, {})
            self.assertIn("SCF not converged", problems["ground_state"])
        with self.assertRaises(ValueError):  # QE input needs the Espresso calculator
            codes.CODES["qe"].plan(bulk("Al"), None, {})

    def test_relative_pseudo_dir_made_absolute(self):
        # pw.x runs inside each run folder: a relative pseudo_dir must not be resolved there
        from ase.calculators.espresso import Espresso, EspressoProfile
        from elastic_ase.rundir import set_directory
        calc = Espresso(profile=EspressoProfile(command="pw.x", pseudo_dir="pseudo"),
                        pseudopotentials={"Al": "Al.UPF"},
                        input_data={"control": {"pseudo_dir": "pp"}})
        set_directory(calc, "runs/ground_state")
        self.assertEqual(calc.profile.pseudo_dir, os.path.abspath("pseudo"))
        self.assertEqual(calc.parameters["input_data"]["control"]["pseudo_dir"], os.path.abspath("pp"))

    def test_relaxation_default_and_off(self):
        # 2-atom hcp cell: the atoms move under basal shear (C66 = (C11 - C12)/2), so
        # relaxed-ion and clamped-ion constants differ
        def c66(relaxation):
            atoms = bulk("Cu", "hcp", a=2.54, c=4.15)
            atoms.calc = EMT()
            ec = ElasticASE(atoms, verbose=False, relaxation=relaxation, fmax=1e-4)
            self.assertEqual(ec.relaxes_atoms, relaxation)
            c = ec.calculate()["standard_constants"]
            return (c["C11"] - c["C12"]) / 2
        relaxed, clamped = c66(True), c66(False)
        self.assertLess(relaxed, clamped - 10.0)  # relaxation lowers the shear stiffness (52 vs 69 GPa)
        al = ElasticASE(bulk("Al", "fcc", a=A0), verbose=False)  # 1 atom: nothing to relax
        self.assertTrue(al.relaxation)
        self.assertFalse(al.relaxes_atoms)
        self.assertEqual(al.optimizer.__name__, "BFGS")

    def test_saved_energies_follow_relaxation_and_structure(self):
        # energies saved with one relaxation setting are not reused for another, nor for a
        # different structure in the same cell
        def hcp():
            atoms = bulk("Cu", "hcp", a=2.54, c=4.15)
            atoms.calc = EMT()
            return atoms
        with tempfile.TemporaryDirectory() as tmp:
            clamped = ElasticASE(hcp(), verbose=False, relaxation=False).calculate(workdir=tmp)
            relaxed = ElasticASE(hcp(), verbose=False).calculate(workdir=tmp)
            in_memory = ElasticASE(hcp(), verbose=False).calculate()
            np.testing.assert_allclose(relaxed.stiffness_matrix, in_memory.stiffness_matrix, atol=1e-6)
            self.assertGreater(abs(relaxed.standard_constants["C11"] - clamped.standard_constants["C11"]), 5)
            ec = ElasticASE(hcp(), verbose=False, relaxation=False)
            self.assertFalse(ec.is_complete(tmp))  # replaced by the relaxed ones: one energy per run
            moved = hcp()
            moved.positions[1] += [0.02, 0.0, 0.0]
            ec = ElasticASE(moved, verbose=False, relaxation=False, reorient=False, symmetry="triclinic")
            self.assertEqual(len(ec.missing(tmp)), len(ec.runs()))

    def test_ground_state_relaxed_too(self):
        # leftover forces in the input must not bias the constants: the unstrained cell is
        # relaxed with the same criterion as the strained ones
        atoms = bulk("Cu", "fcc", a=3.59, cubic=True)
        atoms.calc = EMT()
        reference = ElasticASE(atoms, verbose=False, symmetry="cubic", fmax=1e-4).calculate()
        rattled = atoms.copy()
        rattled.calc = EMT()
        rattled.positions += np.random.default_rng(1).normal(0, 0.01, rattled.positions.shape)
        rattled.set_cell(atoms.cell, scale_atoms=False)
        result = ElasticASE(rattled, verbose=False, symmetry="cubic", reorient=False, fmax=1e-4).calculate()
        for name in ("C11", "C12", "C44"):
            self.assertAlmostEqual(result.standard_constants[name], reference.standard_constants[name], delta=0.5)

    def test_strains_must_differ(self):
        with self.assertRaises(ValueError):
            ElasticASE(bulk("Al", "fcc", a=A0), verbose=False, strains=[-0.01, 0.005, 0.0050000004, 0.01])

    def test_cli_strains_relaxation_and_csv(self):
        from elastic_ase.cli import main
        with tempfile.TemporaryDirectory() as tmp:
            path, csv_path = os.path.join(tmp, "POSCAR"), os.path.join(tmp, "out.csv")
            bulk("Cu", "hcp", a=2.54, c=4.15).write(path, format="vasp")
            self.assertEqual(main([path, "-q", "--no-relaxation", "--strains", "-0.01,0,0.01",
                                   "-o", csv_path]), 0)
            with open(csv_path) as f:
                header, row = [line.strip().split(",") for line in f][:2]
            atoms = bulk("Cu", "hcp", a=2.54, c=4.15)
            atoms.calc = EMT()
            api = ElasticASE(atoms, verbose=False, relaxation=False, strains=[-0.01, 0, 0.01]).calculate()
            self.assertAlmostEqual(float(row[header.index("C11_GPa")]), api.stiffness_matrix[0, 0], delta=0.01)
            self.assertEqual(main([path, "-q", "--strains", "0.01"]), 1)  # clean error, no traceback

    def test_cli_calculator_option(self):
        from elastic_ase.cli import make_calculator
        self.assertIsInstance(make_calculator("emt"), EMT)
        self.assertIsInstance(make_calculator("ase.calculators.emt:EMT"), EMT)
        with self.assertRaises(ValueError):
            make_calculator("lj")
        with self.assertRaises(ImportError):
            make_calculator("no_such_module:Calc")
        # a calculator that needs arguments: clean error and exit code 1, no traceback
        from elastic_ase.cli import main
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "POSCAR")
            bulk("Al", "fcc", a=A0).write(path, format="vasp")
            self.assertEqual(main([path, "-q", "-c", "ase.build:bulk"]), 1)

    def test_deprecated_helpers_match_tables(self):
        from elastic_ase import SYSTEMS, get_deformation_specs, solve_stiffness
        for name, system in SYSTEMS.items():
            self.assertEqual([f for f, _ in get_deformation_specs(name)], [d.folder for d in system.deformations])
            p = np.arange(1, len(system.constants) + 1, dtype=float) * 7.0
            C = system.matrix(p)
            coeffs = {d.folder: 0.5 * np.array(d.strain) @ C @ np.array(d.strain) for d in system.deformations}
            np.testing.assert_allclose(solve_stiffness(name, coeffs), C, atol=1e-9)

if __name__ == "__main__":
    unittest.main()
