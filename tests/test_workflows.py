"""Run folders, saved energies, failure handling and pre/post-processing (no DFT code needed)."""
import json
import os
import tempfile
import unittest

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

from elastic_ase import ElasticASE
from elastic_ase.calculator import _vasp_converged

def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


def _load_json(path):
    with open(path) as f:
        return json.load(f)


A0 = 3.9943  # EMT equilibrium lattice constant of Al
STRAINS = [-0.004, -0.002, 0.0, 0.002, 0.004]


class FileEMT(EMT):
    """EMT that behaves like a file-based DFT calculator: it has a ``directory``, records
    where each run happened, and can fail or report an unconverged SCF for chosen runs."""

    def __init__(self, fail=(), unconverged=()):
        super().__init__()
        self.directory = "."
        self.fail, self.unconverged = set(fail), set(unconverged)
        self.ran_in = []
        self.converged = True

    def calculate(self, atoms=None, properties=("energy",), system_changes=None):
        run = os.path.relpath(self.directory, self.root) if hasattr(self, "root") else self.directory
        self.ran_in.append(run)
        if run in self.fail:
            raise RuntimeError("VASP crashed")
        super().calculate(atoms, properties, system_changes or ["positions", "cell"])
        self.converged = run not in self.unconverged
        with open(os.path.join(self.directory, "OUT"), "w") as f:
            f.write(str(self.results["energy"]))


def al():
    return bulk("Al", "fcc", a=A0, cubic=True)


def reference():
    atoms = al()
    atoms.calc = EMT()
    return ElasticASE(atoms, strains=STRAINS, verbose=False).calculate()["stiffness_matrix"]


class TestWorkflows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.C_ref = reference()

    def test_runs_skip_repeated_zero_strain(self):
        atoms = al()
        atoms.calc = EMT()
        el = ElasticASE(atoms, strains=STRAINS, verbose=False)
        keys = [k for k, _ in el.runs()]
        self.assertEqual(len(keys), 1 + 3 * 4)  # ground state + 3 cubic modes x 4 nonzero strains
        self.assertEqual(keys[0], "ground_state")
        self.assertEqual(len(set(keys)), len(keys))

    def test_workdir_folders_failures_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            fail_key = "C44/strain_+0.002000"
            calc = FileEMT(fail=[fail_key], unconverged=["C11_C12_I/strain_-0.004000"])
            calc.root = tmp
            atoms = al()
            atoms.calc = calc
            el = ElasticASE(atoms, strains=STRAINS, verbose=False)
            with self.assertRaises(RuntimeError) as ctx:
                el.calculate(workdir=tmp)
            message = str(ctx.exception)
            self.assertIn(fail_key, message)
            self.assertIn("C11_C12_I/strain_-0.004000", message)
            self.assertIn("2 of 13", message)
            # every run had its own folder; the calculator's directory is restored afterwards
            self.assertEqual(len(set(calc.ran_in)), 13)
            self.assertTrue(os.path.exists(os.path.join(tmp, "C11_C12_I", "strain_+0.004000", "OUT")))
            self.assertEqual(calc.directory, ".")
            saved = _load_json(os.path.join(tmp, "energies.json"))
            self.assertEqual(len(saved), 11)

            # second call: only the two failed runs are repeated
            calc.fail, calc.unconverged, calc.ran_in = set(), set(), []
            res = ElasticASE(atoms, strains=STRAINS, verbose=False).calculate(workdir=tmp)
            self.assertEqual(sorted(calc.ran_in), sorted([fail_key, "C11_C12_I/strain_-0.004000"]))
            np.testing.assert_allclose(res["stiffness_matrix"], self.C_ref, atol=1e-8)

    def test_saved_energies_of_another_structure_are_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            atoms = al()
            atoms.calc = EMT()
            ElasticASE(atoms, strains=STRAINS, verbose=False).calculate(workdir=tmp)
            other = bulk("Al", "fcc", a=4.05, cubic=True)
            calc = FileEMT()
            calc.root = tmp
            other.calc = calc
            ElasticASE(other, strains=STRAINS, verbose=False).calculate(workdir=tmp)
            self.assertEqual(len(calc.ran_in), 13)  # nothing reused: the cells differ

    def test_write_inputs_then_read_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            atoms = al()  # no calculator needed for pre/post-processing
            el = ElasticASE(atoms, strains=STRAINS, verbose=False)
            folders = el.write_inputs(tmp, format="vasp")
            self.assertEqual(len(folders), 13)
            self.assertTrue(os.path.exists(os.path.join(tmp, "C44", "strain_-0.004000", "POSCAR")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "elastic_ase_runs.json")))

            # "run" every folder outside elastic_ase and leave an output file behind
            for folder in folders:
                structure = read(os.path.join(folder, "POSCAR"))
                structure.calc = EMT()
                energy = structure.get_potential_energy()
                structure.calc = SinglePointCalculator(structure, energy=energy)
                write(os.path.join(folder, "result.extxyz"), structure)

            res = ElasticASE(al(), strains=STRAINS, verbose=False).read_results(tmp, filename="result.extxyz")
            np.testing.assert_allclose(res["stiffness_matrix"], self.C_ref, atol=1e-6)

            # missing and mixed-up outputs are all reported together
            os.remove(os.path.join(folders[3], "result.extxyz"))
            wrong = read(os.path.join(folders[5], "result.extxyz"))
            wrong.set_cell(wrong.cell * 1.01, scale_atoms=True)
            write(os.path.join(folders[5], "result.extxyz"), wrong)
            with self.assertRaises(RuntimeError) as ctx:
                ElasticASE(al(), strains=STRAINS, verbose=False).read_results(tmp, filename="result.extxyz")
            self.assertIn("2 of 13", str(ctx.exception))
            self.assertIn("missing", str(ctx.exception))
            self.assertIn("different cell", str(ctx.exception))

    def test_status_and_run_pending(self):
        with tempfile.TemporaryDirectory() as tmp:
            calc = FileEMT(fail=["C44/strain_-0.002000"])
            calc.root = tmp
            atoms = al()
            atoms.calc = calc
            el = ElasticASE(atoms, strains=STRAINS, verbose=False)
            self.assertEqual(len(el.status(tmp)["pending"]), 13)  # nothing run yet

            # run only two chosen runs
            state = el.run_pending(tmp, runs=["ground_state", "C44/strain_-0.002000"])
            self.assertEqual(state["done"], ["ground_state"])
            self.assertIn("C44/strain_-0.002000", state["failed"])
            self.assertIn("VASP crashed", state["failed"]["C44/strain_-0.002000"])
            self.assertEqual(len(state["pending"]), 11)
            with self.assertRaises(ValueError):
                el.run_pending(tmp, runs=["no_such_run"])

            # everything not done: the failed one still fails, the pending ones succeed
            self.assertFalse(el.is_complete(tmp))
            self.assertEqual(len(el.missing(tmp)), 12)
            state = el.run_pending(tmp)
            self.assertEqual((len(state["done"]), len(state["failed"]), len(state["pending"])), (12, 1, 0))
            self.assertEqual(el.missing(tmp), ["C44/strain_-0.002000"])

            # fixed: re-evaluate only the failed run, then the fit works without new runs
            calc.fail, calc.ran_in = set(), []
            state = el.run_pending(tmp, runs=list(state["failed"]))
            self.assertEqual(calc.ran_in, ["C44/strain_-0.002000"])
            self.assertEqual((len(state["done"]), state["failed"], state["pending"]), (13, {}, []))
            self.assertTrue(el.is_complete(tmp))
            self.assertEqual(el.missing(tmp), [])
            calc.ran_in = []
            res = el.calculate(workdir=tmp)
            self.assertEqual(calc.ran_in, [])
            np.testing.assert_allclose(res["stiffness_matrix"], self.C_ref, atol=1e-8)

    def test_status_from_output_files_and_resubmit(self):
        with tempfile.TemporaryDirectory() as tmp:
            el = ElasticASE(al(), strains=STRAINS, verbose=False)
            folders = el.write_inputs(tmp)
            for folder in folders[:10]:  # three runs not done yet
                structure = read(os.path.join(folder, "POSCAR"))
                structure.calc = EMT()
                structure.calc = SinglePointCalculator(structure, energy=structure.get_potential_energy())
                write(os.path.join(folder, "out.extxyz"), structure)
            _write(os.path.join(folders[10], "out.extxyz"), "garbage")  # a crashed run
            state = el.status(tmp, filename="out.extxyz")
            self.assertEqual((len(state["done"]), len(state["failed"]), len(state["pending"])), (10, 1, 2))
            redo = el.missing(tmp, filename="out.extxyz")
            self.assertEqual(sorted(redo), sorted(list(state["failed"]) + state["pending"]))
            self.assertFalse(el.is_complete(tmp, filename="out.extxyz"))
            rewritten = el.write_inputs(tmp, runs=redo)
            self.assertEqual(sorted(rewritten), sorted(os.path.join(tmp, k) for k in redo))

    def test_vasp_input_copies_user_files_and_checks_potcar(self):
        from ase import Atoms
        with tempfile.TemporaryDirectory() as tmp:
            incar, kpoints, potcar = (os.path.join(tmp, n) for n in ("INCAR", "KPOINTS", "POTCAR"))
            _write(incar, "ENCUT = 520\nISIF = 2\nIBRION = -1\n")
            _write(kpoints, "auto\n0\nGamma\n12 12 12\n")
            _write(potcar, "  PAW_PBE Al 04Jan2001\n   TITEL  = PAW_PBE Al 04Jan2001\n")
            ec = ElasticASE(al(), strains=STRAINS, verbose=False)
            runs = os.path.join(tmp, "runs")
            folders = ec.vasp_input(runs, incar=incar, kpoints=kpoints, potcar=potcar)
            self.assertEqual(len(folders), 13)
            for folder in folders:
                self.assertEqual(sorted(os.listdir(folder)), ["INCAR", "KPOINTS", "POSCAR", "POTCAR"])
            with open(os.path.join(folders[5], "INCAR")) as copied, open(incar) as original:
                self.assertEqual(copied.read(), original.read())
            strained = read(os.path.join(runs, "C44", "strain_+0.004000", "POSCAR"))
            self.assertAlmostEqual(strained.cell[1, 2], A0 * 0.002, places=8)  # F = I + eps, eps_yz = 0.004/2

            # resubmit only chosen runs
            again = ec.vasp_input(runs, incar=incar, kpoints=kpoints, potcar=potcar, runs=["C44/strain_-0.004000"])
            self.assertEqual(again, [os.path.join(runs, "C44/strain_-0.004000")])

            # POTCAR for another element, a missing file, atoms not grouped by element
            _write(potcar, "   TITEL  = PAW_PBE Cu 22Jun2005\n")
            with self.assertRaises(ValueError):
                ec.vasp_input(runs, potcar=potcar)
            with self.assertRaises(FileNotFoundError):
                ec.vasp_input(runs, incar=os.path.join(tmp, "nope"))
            mixed = Atoms("CuAuCu", scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0, 0]],
                          cell=[3, 3, 6], pbc=True)
            with self.assertRaises(ValueError):
                ElasticASE(mixed, symmetry="triclinic", reorient=False, strains=STRAINS,
                           verbose=False).vasp_input(os.path.join(tmp, "mixed"))

    def test_vasp_status_and_output_are_the_outcar_versions(self):
        with tempfile.TemporaryDirectory() as tmp:
            ec = ElasticASE(al(), strains=STRAINS, verbose=False)
            ec.vasp_input(tmp)
            self.assertEqual(len(ec.vasp_status(tmp)["pending"]), 13)  # no OUTCARs yet
            self.assertFalse(ec.is_complete(tmp, filename="OUTCAR"))
            with self.assertRaises(RuntimeError):
                ec.vasp_output(tmp)

    def test_code_setting_and_run_jobs(self):
        import sys
        with tempfile.TemporaryDirectory() as tmp:
            program = os.path.join(tmp, "fake_dft.py")  # stands in for vasp_std / pw.x
            _write(program, (
                "import os, sys\n"
                "from ase.io import read, write\n"
                "from ase.calculators.emt import EMT\n"
                "from ase.calculators.singlepoint import SinglePointCalculator\n"
                "if os.path.exists('CRASH'):\n"
                "    sys.exit('segmentation fault')\n"
                "a = read('POSCAR')\n"
                "a.calc = EMT()\n"
                "e = a.get_potential_energy()\n"
                "a.calc = SinglePointCalculator(a, energy=e)\n"
                "write('out.extxyz', a)\n"))
            runs = os.path.join(tmp, "runs")
            command = f'"{sys.executable}" "{program}"'
            ec = ElasticASE(al(), strains=STRAINS, verbose=False, output="out.extxyz", command=command)
            folders = ec.write_inputs(runs, files={"INCAR": program})  # any extra file is copied
            self.assertTrue(os.path.exists(os.path.join(folders[0], "INCAR")))
            crash = os.path.join(runs, "C44", "strain_+0.004000", "CRASH")
            _write(crash, "")

            state = ec.run_jobs(runs)
            self.assertEqual((len(state["done"]), len(state["pending"])), (12, 1))
            self.assertNotEqual(ec.exit_codes["C44/strain_+0.004000"], 0)
            self.assertFalse(ec.is_complete(runs))
            with open(os.path.join(runs, "C44", "strain_+0.004000", "elastic_ase_run.log")) as f:
                self.assertIn("segmentation fault", f.read())

            os.remove(crash)
            ec.run_jobs(runs)
            self.assertEqual(list(ec.exit_codes), ["C44/strain_+0.004000"])  # only the missing one
            self.assertTrue(ec.is_complete(runs))
            res = ec.read_results(runs)
            np.testing.assert_allclose(res["stiffness_matrix"], self.C_ref, atol=1e-6)

            with self.assertRaises(ValueError):  # no command and no code
                ElasticASE(al(), strains=STRAINS, verbose=False, output="x").run_jobs(runs)

    def test_code_names_and_defaults(self):
        ec = ElasticASE(al(), strains=STRAINS, verbose=False, code="VASP")
        self.assertEqual((ec.code, ec._output_name()), ("vasp", "OUTCAR"))
        self.assertEqual(ElasticASE(al(), verbose=False, code="espresso")._output_name(), "espresso.pwo")
        self.assertEqual(ElasticASE(al(), verbose=False, code="vasp", output="vasprun.xml")._output_name(), "vasprun.xml")
        self.assertEqual(ec._output_name(code="qe"), "espresso.pwo")
        with self.assertRaises(ValueError):
            ElasticASE(al(), verbose=False, code="castepp")
        with tempfile.TemporaryDirectory() as tmp:
            ec.write_inputs(tmp)  # code="vasp": POSCAR folders
            self.assertTrue(os.path.exists(os.path.join(tmp, "C44", "strain_+0.002000", "POSCAR")))
            self.assertEqual(len(ec.status(tmp)["pending"]), 13)  # reads OUTCARs (none yet)

    def test_write_inputs_uses_calculator_writer(self):
        class Writer(EMT):
            directory = "."

            def write_input(self, atoms, properties=None, system_changes=None):
                write(os.path.join(self.directory, "INPUT"), atoms, format="vasp")

        with tempfile.TemporaryDirectory() as tmp:
            atoms = al()
            atoms.calc = Writer()
            folders = ElasticASE(atoms, strains=STRAINS, verbose=False).write_inputs(tmp)
            self.assertTrue(all(os.path.exists(os.path.join(f, "INPUT")) for f in folders))
            self.assertEqual(atoms.calc.directory, ".")

    def test_vasp_convergence_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            good, bad, other = (os.path.join(tmp, n) for n in ("OUTCAR", "OUTCAR_bad", "vasprun.xml"))
            _write(good, "------- aborting loop because EDIFF is reached ----\n")
            _write(bad, "------- aborting loop EDIFF was not reached (unconverged) ----\n")
            _write(other, "<xml/>")
            self.assertTrue(_vasp_converged(good))
            self.assertFalse(_vasp_converged(bad))
            self.assertIsNone(_vasp_converged(other))

    def test_calculate_needs_a_calculator(self):
        with self.assertRaises(ValueError):
            ElasticASE(al(), strains=STRAINS, verbose=False).calculate()


if __name__ == "__main__":
    unittest.main()
