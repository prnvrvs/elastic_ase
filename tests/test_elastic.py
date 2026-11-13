import os
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

        # 1. Rhombohedral cell check
        bi = bulk("Bi", "rhombohedral", a=4.75, alpha=57.2)
        warn_rhombo = validate_orientation(bi, "trigonal_2", verbose=False)
        self.assertTrue(any("Rhombohedral cell detected" in w for w in warn_rhombo))

        # 2. Hexagonal tilted c-axis check (c along x)
        co = bulk("Co", "hcp", a=2.5, c=4.0)
        c_orig = co.get_cell()
        co.set_cell([c_orig[2], c_orig[1], c_orig[0]], scale_atoms=True)
        warn_hex = validate_orientation(co, "hexagonal", verbose=False)
        self.assertTrue(any("c-axis is tilted" in w for w in warn_hex))

        # 3. Monoclinic with unique axis c (gamma != 90)
        mono_c = Atoms("C", cell=[[3.0, 0, 0], [1.0, 4.0, 0], [0, 0, 5.0]], pbc=True)
        warn_mono = validate_orientation(mono_c, "monoclinic", verbose=False)
        self.assertTrue(any("unique axis c" in w for w in warn_mono))

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

        # Check directional properties calculated without plotting
        dp = res["directional_properties"]
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

        # Check Universal Anisotropy in averages
        self.assertIn("universal_anisotropy", res["polycrystalline_averages"])
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

if __name__ == "__main__":
    unittest.main()
