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
