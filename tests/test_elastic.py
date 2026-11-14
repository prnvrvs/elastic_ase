import os
import importlib.util
import unittest
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from elastic_ase import ElasticASE

class TestElasticASE(unittest.TestCase):
    def test_cubic_emt(self):
        al = bulk("Al", "fcc", a=4.05, cubic=True)
        al.calc = EMT()

        el = ElasticASE(al, symmetry="auto", verbose=False)
        res = el.calculate()

        self.assertEqual(el.symmetry, "cubic")
        self.assertTrue(el.is_stable)

        C = res["stiffness_matrix"]
        self.assertEqual(C.shape, (6, 6))
        self.assertGreater(C[0, 0], 0)  # C11
        self.assertGreater(C[0, 1], 0)  # C12
        self.assertGreater(C[3, 3], 0)  # C44

        # Cubic symmetry equalities
        np.testing.assert_allclose(C[0, 0], C[1, 1], rtol=1e-3)
        np.testing.assert_allclose(C[0, 0], C[2, 2], rtol=1e-3)
        np.testing.assert_allclose(C[0, 1], C[0, 2], rtol=1e-3)
        np.testing.assert_allclose(C[3, 3], C[4, 4], rtol=1e-3)
        np.testing.assert_allclose(C[3, 3], C[5, 5], rtol=1e-3)

    def test_hexagonal_emt(self):
        # Ni in HCP structure with EMT calculator
        ni_hcp = bulk("Ni", "hcp", a=2.50, c=4.08)
        ni_hcp.calc = EMT()

        el = ElasticASE(ni_hcp, symmetry="auto", verbose=False)
        res = el.calculate()

        self.assertEqual(el.symmetry, "hexagonal")
        self.assertTrue(el.is_stable)

        C = res["stiffness_matrix"]
        self.assertGreater(C[0, 0], 0)  # C11
        self.assertGreater(C[0, 1], 0)  # C12
        self.assertGreater(C[0, 2], 0)  # C13
        self.assertGreater(C[2, 2], 0)  # C33
        self.assertGreater(C[3, 3], 0)  # C44

        # Hexagonal relation: C66 = (C11 - C12) / 2
        np.testing.assert_allclose(C[5, 5], (C[0, 0] - C[0, 1]) / 2.0, rtol=1e-3)

    def test_polycrystalline_averages(self):
        al = bulk("Al", "fcc", a=4.05, cubic=True)
        al.calc = EMT()

        el = ElasticASE(al, symmetry="cubic", verbose=False)
        res = el.calculate()
        avg = res["polycrystalline_averages"]

        # Voigt-Reuss bounds: K_reuss <= K_VRH <= K_voigt
        self.assertLessEqual(avg["K_reuss"], avg["K_VRH"] + 1e-4)
        self.assertLessEqual(avg["K_VRH"], avg["K_voigt"] + 1e-4)
        self.assertGreater(avg["Pugh_ratio"], 0)

    def test_monoclinic_deformation_specs(self):
        from elastic_ase.deformations import get_deformation_specs
        specs = dict(get_deformation_specs("monoclinic"))
        def_func_II = specs["C44_C66_C46_II"]
        # Must strain epsilon_4 and epsilon_6 (indices 3 and 5), NOT epsilon_5 (index 4)
        d1, d2, d3, d4, d5, d6 = def_func_II(0.01)
        self.assertEqual(d4, 0.01)
        self.assertEqual(d5, 0.0)
        self.assertEqual(d6, -0.01)

    def test_monoclinic_solver_signs(self):
        from elastic_ase.solvers import solve_stiffness
        # Create coefficients where C15 > 0 and C46 > 0
        # H = 0.5*(C11 + C55) + C15, I = 0.5*(C11 + C55) - C15 -> H - I = 2*C15
        coeffs = {
            "C11_C22_C12_I": 50.0, "C11_C22_C12_II": 150.0,
            "C22_C33_C23_I": 50.0, "C22_C33_C23_II": 150.0,
            "C11_C33_C13_I": 50.0, "C11_C33_C13_II": 150.0,
            "C44": 25.0,
            "C11_C55_C15_I": 120.0,  # H
            "C11_C55_C15_II": 80.0,  # I -> C15 = (120 - 80)/2 = +20
            "C22_C55_C25": 100.0,
            "C33_C55_C35": 100.0,
            "C44_C66_C46_I": 90.0,   # L
            "C44_C66_C46_II": 70.0,  # M -> C46 = (90 - 70)/2 = +10
        }
        C = solve_stiffness("monoclinic", coeffs)
        self.assertAlmostEqual(C[0, 4], 20.0)  # C15
        self.assertAlmostEqual(C[3, 5], 10.0)  # C46

    def test_trigonal_2_c15_sign(self):
        from elastic_ase.solvers import solve_stiffness
        # G = 0.5*(C11 + C44) - C15 -> C15 = 0.5*(C11 + C44) - G
        # If C11 = 100, C44 = 40, and C15 = +10 -> G = 70 - 10 = 60
        coeffs = {
            "C11_C12_I": 120.0, "C11_C12_II": 80.0,  # C11 = 100, C12 = 20
            "C11_C33_C13_I": 110.0, "C11_C33_C13_II": 90.0,
            "C11_C44_C14_I": 75.0, "C11_C44_C14_II": 65.0,  # C44 = 40, C14 = 5
            "C11_C44_C15": 60.0  # G -> C15 = 0.5*(100 + 40) - 60 = +10
        }
        C = solve_stiffness("trigonal_2", coeffs)
        self.assertAlmostEqual(C[0, 4], 10.0)   # C15
        self.assertAlmostEqual(C[1, 4], -10.0)  # C25 = -C15
        self.assertAlmostEqual(C[3, 5], -10.0)  # C46 = -C15

    def test_orientation_warnings(self):
        from elastic_ase.calculator import validate_orientation
        from ase import Atoms
        import warnings

        def check(atoms, symmetry):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                return validate_orientation(atoms, symmetry, verbose=False)

        # Correctly oriented crystals: no warning
        self.assertEqual(check(bulk("Cu", "fcc", a=3.6, cubic=True), "cubic"), [])
        self.assertEqual(check(bulk("Cu", "fcc", a=3.6), "cubic"), [])  # primitive cell, cubic axes along x, y, z
        self.assertEqual(check(bulk("Co", "hcp", a=2.5, c=4.0), "hexagonal"), [])

        # 1. Rhombohedral Bi as built by ASE: 3-fold axis not along z
        bi = bulk("Bi", "rhombohedral", a=4.75, alpha=57.2)
        self.assertEqual(len(check(bi, "trigonal_2")), 1)

        # 2. Hexagonal with the c axis along x
        co = bulk("Co", "hcp", a=2.5, c=4.0)
        c_orig = co.get_cell()
        co.set_cell([c_orig[2], c_orig[1], c_orig[0]], scale_atoms=True)
        self.assertEqual(len(check(co, "hexagonal")), 1)

        # 3. Monoclinic with unique axis c (gamma != 90): the table expects the unique axis along y
        mono_c = Atoms("C", cell=[[3.0, 0, 0], [1.0, 4.0, 0], [0, 0, 5.0]], pbc=True)
        self.assertEqual(len(check(mono_c, "monoclinic")), 1)

        # 4. Cubic crystal rotated about z: cubic axes no longer along x, y
        cu = bulk("Cu", "fcc", a=3.6, cubic=True)
        cu.rotate(30, "z", rotate_cell=True)
        self.assertEqual(len(check(cu, "cubic")), 1)
        self.assertEqual(check(cu, "triclinic"), [])  # the full table fits any orientation

    def test_trigonal_auto_routes_to_trigonal_2(self):
        from elastic_ase.calculator import auto_detect_symmetry
        bi = bulk("Bi", "rhombohedral", a=4.75, alpha=57.2)
        sym = auto_detect_symmetry(bi)
        self.assertEqual(sym, "trigonal_2")


    def test_poscar_auto_reorientation(self):
        import tempfile
        from ase.io import read as ase_read

        # Standard reference untilted Ni
        ni_ref = bulk("Ni", "hcp", a=2.50, c=4.08)
        ni_ref.calc = EMT()
        el_ref = ElasticASE(ni_ref, symmetry="hexagonal", verbose=False)
        res_ref = el_ref.calculate()

        # Tilted Ni (c axis along x) written to POSCAR format
        ni_tilted = bulk("Ni", "hcp", a=2.50, c=4.08)
        ni_tilted.rotate(90, 'y', rotate_cell=True)

        with tempfile.NamedTemporaryFile(suffix=".vasp", delete=False) as f:
            poscar_path = f.name

        try:
            ni_tilted.write(poscar_path, format="vasp")
            # Read POSCAR
            atoms_from_poscar = ase_read(poscar_path, format="vasp")
            atoms_from_poscar.calc = EMT()

            # ElasticASE should detect orientation, decide reorientation is required, and reorient
            el_poscar = ElasticASE(atoms_from_poscar, verbose=False)
            self.assertTrue(el_poscar.reoriented)
            self.assertEqual(el_poscar.symmetry, "hexagonal")

            res_poscar = el_poscar.calculate()
            # Verified that elastic constants match reference down to numerical precision
            np.testing.assert_allclose(
                res_poscar["stiffness_matrix"],
                res_ref["stiffness_matrix"],
                atol=1e-2
            )
        finally:
            if os.path.exists(poscar_path):
                os.remove(poscar_path)

    def test_check_symmetry(self):
        from elastic_ase import check_symmetry
        al = bulk("Al", "fcc", a=4.05, cubic=True)
        report = check_symmetry(al, verbose=False)
        self.assertEqual(report["crystal_system"], "cubic")
        self.assertEqual(report["recommended"], "cubic")
        self.assertEqual(report["n_constants"], 3)

    def test_tetragonal_solvers(self):
        from elastic_ase.solvers import solve_stiffness
        # Tetragonal 1: 6 constants
        coeffs_t1 = {
            "C11_C12_I": 120.0, "C11_C12_II": 80.0,  # C11 = 100, C12 = 20
            "C11_C33_C13_I": 90.0,
            "C11_C12_C13_C33": 250.0,
            "C44": 35.0,
            "C66": 25.0
        }
        C1 = solve_stiffness("tetragonal_1", coeffs_t1)
        self.assertEqual(C1.shape, (6, 6))
        self.assertAlmostEqual(C1[0, 0], 100.0)  # C11
        self.assertAlmostEqual(C1[0, 1], 20.0)   # C12
        self.assertAlmostEqual(C1[3, 3], 35.0)   # C44
        self.assertAlmostEqual(C1[5, 5], 50.0)   # C66 = 2 * 25

        # Tetragonal 2: adds C16
        coeffs_t2 = dict(coeffs_t1)
        # G = 0.5*C11 + 0.5*C66 + C16 = 50 + 25 + 15 = 90
        coeffs_t2["C11_C66_C16"] = 90.0
        C2 = solve_stiffness("tetragonal_2", coeffs_t2)
        self.assertAlmostEqual(C2[0, 5], 15.0)   # C16
        self.assertAlmostEqual(C2[1, 5], -15.0)  # C26 = -C16
        self.assertAlmostEqual(C2[5, 0], 15.0)
        self.assertAlmostEqual(C2[5, 1], -15.0)

    @unittest.skipUnless(importlib.util.find_spec("elastic_vis"), "elastic-vis not installed")
    def test_elastic_vis_visualization(self):
        import tempfile
        al = bulk("Al", "fcc", a=4.05, cubic=True)
        al.calc = EMT()

        el = ElasticASE(al, symmetry="cubic", verbose=False)
        el.calculate()

        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "test_al")
            # Generate 2D and 3D plots for Young's modulus
            plots = el.visualize(
                properties=["Young"],
                dimensions=["2D", "3D"],
                prefix=prefix,
                output_format="png"
            )
            self.assertEqual(len(plots), 2)
            for p in plots:
                self.assertTrue(os.path.exists(p))
                self.assertGreater(os.path.getsize(p), 0)

            # Test export_matrix
            txt_path = os.path.join(tmpdir, "C.txt")
            el.export_matrix(txt_path)
            self.assertTrue(os.path.exists(txt_path))
            loaded = np.loadtxt(txt_path)
            np.testing.assert_allclose(loaded, el.stiffness)

    def test_directional_properties_and_compliance(self):
        al = bulk("Al", "fcc", a=4.05, cubic=True)
        al.calc = EMT()

        el = ElasticASE(al, symmetry="cubic", verbose=False)
        res = el.calculate()

        # Check compliance matrix S = C^-1
        S = el.compliance
        np.testing.assert_allclose(S @ el.stiffness, np.eye(6), atol=1e-5)

        # Check Universal Anisotropy in averages
        self.assertIn("universal_anisotropy", res["polycrystalline_averages"])

        # Directional properties need the optional elastic-vis; without it they are None
        dp = res["directional_properties"]
        if importlib.util.find_spec("elastic_vis") is None:
            self.assertIsNone(dp)
            return
        self.assertIsNotNone(dp)
        self.assertGreater(dp["E_min"], 0)
        self.assertGreaterEqual(dp["E_max"], dp["E_min"])
        self.assertGreaterEqual(dp["E_anisotropy"], 1.0)
        self.assertGreater(dp["G_min"], 0)
        self.assertGreaterEqual(dp["G_max"], dp["G_min"])
        self.assertIsInstance(dp["is_auxetic"], bool)
        self.assertIsInstance(dp["is_negative_LC"], bool)

        # Check direct directional methods
        # In [001] direction: theta = 0
        E_001 = el.Young(0.0, 0.0)
        self.assertAlmostEqual(E_001, 1.0 / S[2, 2], places=3)
        self.assertGreaterEqual(res["polycrystalline_averages"]["universal_anisotropy"], 0.0)

    def test_cli_execution(self):
        import tempfile
        from elastic_ase.cli import main as cli_main

        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as f:
            cif_path = f.name
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            json_path = f.name

        try:
            bulk("Al", "fcc", a=4.05, cubic=True).write(cif_path)
            ret = cli_main([cif_path, "-c", "emt", "-o", json_path, "-q"])
            self.assertEqual(ret, 0)
            self.assertTrue(os.path.exists(json_path))
            self.assertGreater(os.path.getsize(json_path), 0)
        finally:
            if os.path.exists(cif_path):
                os.remove(cif_path)
            if os.path.exists(json_path):
                os.remove(json_path)

    def test_orientation_aware_auto_mode_arbitrary_tilt(self):
        # 1. Untilted FCC Al reference
        al_ref = bulk("Al", "fcc", a=4.05, cubic=True)
        al_ref.calc = EMT()
        el_ref = ElasticASE(al_ref, symmetry="cubic", verbose=False)
        res_ref = el_ref.calculate()
        C_ref = res_ref["stiffness_matrix"]

        # 2. Arbitrarily tilted FCC Al (rotated by 37 degrees around [1, 2, 3])
        al_tilted = bulk("Al", "fcc", a=4.05, cubic=True)
        al_tilted.rotate(37, [1, 2, 3], rotate_cell=True)
        al_tilted.calc = EMT()

        # Run with reorient=False (orientation-aware auto mode) and frame="standard"
        el_tilted = ElasticASE(al_tilted, symmetry="auto", reorient=False, frame="standard", verbose=False)
        self.assertTrue(el_tilted.is_auto_oriented)
        res_tilted = el_tilted.calculate()

        C_std = res_tilted["stiffness_standard"]
        # The transformed standard stiffness should match reference C_ref within 1.5 GPa
        np.testing.assert_allclose(C_std[0, 0], C_ref[0, 0], atol=1.5)
        np.testing.assert_allclose(C_std[0, 1], C_ref[0, 1], atol=1.5)
        np.testing.assert_allclose(C_std[3, 3], C_ref[3, 3], atol=1.5)

        # Polycrystalline VRH moduli must be rotationally invariant
        np.testing.assert_allclose(
            res_tilted["polycrystalline_averages"]["K_VRH"],
            res_ref["polycrystalline_averages"]["K_VRH"],
            atol=0.5
        )

    def test_standard_frame_after_reorientation(self):
        # Default reorient=True: a rotated cubic cell is standardized first, so the standard
        # frame must equal the unrotated result (the rotation must not be applied twice).
        a0 = 3.9943  # EMT equilibrium lattice constant of Al (stress-free reference)
        strains = [-0.002, -0.001, 0.0, 0.001, 0.002]
        al_ref = bulk("Al", "fcc", a=a0, cubic=True)
        al_ref.calc = EMT()
        C_ref = ElasticASE(al_ref, symmetry="cubic", strains=strains, verbose=False).calculate()["stiffness_matrix"]

        al_rot = bulk("Al", "fcc", a=a0, cubic=True)
        al_rot.rotate(37, [1, 2, 3], rotate_cell=True)
        al_rot.calc = EMT()
        el = ElasticASE(al_rot, frame="standard", strains=strains, verbose=False)
        self.assertTrue(el.reoriented)
        res = el.calculate()
        np.testing.assert_allclose(res["stiffness_standard"], C_ref, atol=0.05)
        np.testing.assert_allclose(res["stiffness_matrix"], C_ref, atol=0.05)
        for name in ("C11", "C12", "C44"):
            self.assertAlmostEqual(res["standard_constants"][name], C_ref[{"C11": 0, "C12": 0, "C44": 3}[name],
                                   {"C11": 0, "C12": 1, "C44": 3}[name]], delta=0.05)

    def test_reorientation_keeps_cell_and_atom_data(self):
        # Reorientation is a rigid rotation: the user's cell (primitive or supercell), the atoms
        # and per-atom data stay as given, so calculator settings (k-points, magmoms) remain valid.
        from ase.constraints import FixAtoms
        cu = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((2, 1, 1))
        cu.rotate(30, "z", rotate_cell=True)
        cu.set_initial_magnetic_moments([0.5] * len(cu))
        cu.set_tags(range(len(cu)))
        cu.set_constraint(FixAtoms([0]))
        cu.calc = EMT()
        el = ElasticASE(cu, verbose=False)
        self.assertTrue(el.reoriented)
        self.assertEqual(len(el.atoms), len(cu))
        np.testing.assert_allclose(el.atoms.cell.cellpar(), cu.cell.cellpar(), atol=1e-8)
        np.testing.assert_allclose(el.atoms.get_initial_magnetic_moments(), 0.5)
        np.testing.assert_array_equal(el.atoms.get_tags(), cu.get_tags())
        self.assertEqual(len(el.atoms.constraints), 1)
        np.testing.assert_allclose(np.abs(el.atoms.cell[0]), [7.2, 0, 0], atol=1e-8)  # cubic axes along x, y, z

        prim = bulk("Cu", "fcc", a=3.6)
        prim.calc = EMT()
        self.assertEqual(len(ElasticASE(prim, verbose=False).atoms), 1)

    def test_stability_independent_of_orientation(self):
        from elastic_ase import check_mechanical_stability
        from elastic_ase.auto import bond_matrix
        C = np.zeros((6, 6))
        C[:3, :3] = 60.0
        np.fill_diagonal(C[:3, :3], 100.0)
        C[3, 3] = C[4, 4] = C[5, 5] = -5.0  # unstable: C44 < 0
        rng = np.random.default_rng(0)
        for _ in range(50):
            q, _r = np.linalg.qr(rng.normal(size=(3, 3)))
            q *= np.sign(np.linalg.det(q))
            m = bond_matrix(q)
            self.assertFalse(check_mechanical_stability(m @ C @ m.T, "cubic")[0])

    def test_needs_three_strains(self):
        al = bulk("Al", "fcc", a=4.05, cubic=True)
        al.calc = EMT()
        with self.assertRaises(ValueError):
            ElasticASE(al, strains=[0.01], verbose=False)
        with self.assertRaises(ValueError):
            ElasticASE(al, strains=[0.01, 0.01, -0.01], verbose=False)

    def test_laue_classes_of_all_space_groups(self):
        from elastic_ase.symmetry import classify
        expected = {1: "-1", 2: "-1", 3: "2/m", 15: "2/m", 16: "mmm", 74: "mmm", 75: "4/m", 88: "4/m",
                    89: "4/mmm", 142: "4/mmm", 143: "-3", 148: "-3", 149: "-3m", 167: "-3m",
                    168: "6/m", 176: "6/m", 177: "6/mmm", 194: "6/mmm", 195: "m-3", 205: "m-3",
                    206: "m-3", 207: "m-3m", 225: "m-3m", 230: "m-3m"}
        for number, laue in expected.items():
            self.assertEqual(classify(number)[1], laue, number)
        self.assertEqual(classify(176)[2], "hexagonal")   # 6/m and 6/mmm: same 5-constant table
        self.assertEqual(classify(205)[2], "cubic")       # m-3 and m-3m: same 3-constant table

    def test_symmetry_classification_and_report(self):
        from elastic_ase.symmetry import classify, check_symmetry, TRIGONAL_31M

        # Test space group classification
        sys_c, laue_c, name_c = classify(225)
        self.assertEqual((sys_c, laue_c, name_c), ("cubic", "m-3m", "cubic"))

        sys_h, laue_h, name_h = classify(194)
        self.assertEqual((sys_h, laue_h, name_h), ("hexagonal", "6/mmm", "hexagonal"))

        # Test trigonal_31m routing to trigonal_2
        for sg in TRIGONAL_31M:
            _, _, named = classify(sg)
            self.assertEqual(named, "trigonal_2")

        # Test SymmetryReport summary method
        al = bulk("Al", "fcc", a=4.05, cubic=True)
        report = check_symmetry(al, verbose=False)
        self.assertEqual(report.crystal_system, "cubic")
        self.assertEqual(report.laue_class, "m-3m")
        self.assertEqual(report.recommended, "cubic")
        self.assertTrue(report.recommended_works)
        summary_str = report.summary()
        self.assertIn("Fm-3m", summary_str)
        self.assertIn("3 independent", summary_str)

    def test_cli_frame_and_check(self):
        import tempfile
        import json
        from elastic_ase.cli import main as cli_main

        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as f:
            cif_path = f.name
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            json_path = f.name

        try:
            bulk("Al", "fcc", a=4.05, cubic=True).write(cif_path)
            # Test --check mode
            ret_check = cli_main([cif_path, "--check"])
            self.assertEqual(ret_check, 0)

            # Test --frame standard mode
            ret_calc = cli_main([cif_path, "-c", "emt", "--frame", "standard", "-o", json_path, "-q"])
            self.assertEqual(ret_calc, 0)
            with open(json_path) as jf:
                data = json.load(jf)
            self.assertEqual(data["frame"], "standard")
            self.assertIn("stiffness_standard", data)
            self.assertEqual(len(data["stiffness_standard"]), 6)
        finally:
            if os.path.exists(cif_path):
                os.remove(cif_path)
            if os.path.exists(json_path):
                os.remove(json_path)

if __name__ == "__main__":
    unittest.main()
