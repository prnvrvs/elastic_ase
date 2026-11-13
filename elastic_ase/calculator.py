"""The main class, ElasticConstants (= ElasticASE).

It sets up the symmetry and the runs, gets the energies (attached calculator or run folders) and
hands them to :mod:`elastic_ase.analysis`. Run folders are handled by :mod:`elastic_ase.rundir`,
DFT-code specifics by :mod:`elastic_ase.codes`, the report and plots by
:mod:`elastic_ase.report` and :mod:`elastic_ase.plotting`.
"""
import json
import logging
import os
import shutil
import subprocess
import warnings

import numpy as np
from ase import Atoms
from ase.io import write as ase_write

from elastic_ase import analysis, auto, codes, orientation, plotting, report, rundir
from elastic_ase._spglib import call as spg_call
from elastic_ase.deformations import build_strain_matrix
from elastic_ase.rundir import GROUND_STATE, strained_key
from elastic_ase.systems import SYSTEMS, _ALIASES as SYSTEM_ALIASES

# Names kept importable from here for compatibility
from elastic_ase.analysis import EV_TO_GPA, RESIDUAL_STRESS_TOL  # noqa: F401
from elastic_ase.codes import CODES, qe_converged as _qe_converged, vasp_converged as _vasp_converged  # noqa: F401
from elastic_ase.rundir import ENERGY_FILE, FAILED_FILE  # noqa: F401

logger = logging.getLogger("elastic_ase")
DEFAULT_STRAINS = (-0.01, -0.005, 0.0, 0.005, 0.01)


def auto_detect_symmetry(atoms: Atoms, symprec: float = 1e-3) -> str:
    """
    Auto-detect crystal symmetry family using spglib.
    Returns named workflow string ('cubic', 'hexagonal', 'trigonal_2', etc.).
    Falls back to 'triclinic' (all 21 constants, valid for any crystal) if spglib fails.
    """
    cell_data = (np.asarray(atoms.get_cell(), float), atoms.get_scaled_positions(), atoms.numbers)
    dataset = spg_call("get_symmetry_dataset", cell_data, symprec=symprec)
    if dataset is None:
        warnings.warn(
            f"spglib could not determine the symmetry (symprec={symprec}); using the triclinic "
            "workflow (21 constants), which is valid for any crystal."
        )
        return "triclinic"

    spg = dataset.number if hasattr(dataset, "number") else dataset["number"]
    if 195 <= spg <= 230:
        return "cubic"
    elif 168 <= spg <= 194:
        return "hexagonal"
    elif 143 <= spg <= 167:
        # Route all trigonal groups to trigonal_2 because it simultaneously fits C14 and C15,
        # avoiding silent errors in '31m'-type groups and for arbitrary in-plane rotation about z.
        return "trigonal_2"
    elif 89 <= spg <= 142:
        return "tetragonal_1"
    elif 75 <= spg <= 88:
        return "tetragonal_2"
    elif 16 <= spg <= 74:
        return "orthogonal"
    elif 3 <= spg <= 15:
        return "monoclinic"
    else:
        return "triclinic"


def validate_orientation(atoms: Atoms, symmetry: str, verbose: bool = True, symprec: float = 1e-3):
    """
    Check that the named workflow can measure this crystal in its current orientation.

    The named tables assume the crystal axes along fixed Cartesian directions (e.g. cubic axes
    along x, y, z; the 3-, 4- or 6-fold axis along z; the monoclinic unique axis along y). The check
    is exact: the table must be able to represent every stiffness matrix the crystal's point
    group allows in this orientation. Returns a list of warning messages (empty if it fits).
    """
    sym = symmetry.lower().strip()
    if sym not in SYSTEMS:
        return []
    try:
        basis = auto.invariant_basis(auto.crystal_rotations(atoms, symprec=symprec))
    except Exception:
        return []
    if orientation.table_covers(SYSTEMS[sym], basis):
        return []

    message = (
        f"The '{sym}' deformations cannot measure this crystal in its current orientation "
        f"(its symmetry allows {len(basis)} independent constants that the '{sym}' table does not "
        "cover), so the constants would be wrong. Use the default reorient=True, or symmetry='auto' "
        "(works in any orientation), or rotate the structure into the standard setting."
    )
    warnings.warn(message)
    if verbose:
        print(f"[ElasticASE Warning] {message}")
    return [message]


def standardize_and_reorient(atoms: Atoms, symprec: float = 1e-3, verbose: bool = True):
    """
    Rigidly rotate the structure so that its crystal axes are in the standard setting
    (cubic axes along x, y, z; the 3-, 4- or 6-fold axis along z; the monoclinic unique axis b
    along y).

    Only the orientation changes: the cell (primitive, conventional or supercell), the number
    and order of atoms and all per-atom data (magnetic moments, charges, tags, constraints) are
    kept, so calculator settings such as the k-point mesh remain valid.
    Returns (atoms, rotated).
    """
    try:
        rotation, number = auto.standard_setting(atoms, symprec=symprec)
    except Exception as e:
        if verbose:
            print(f"[ElasticASE] Notice: reorientation skipped ({e})")
        return atoms, False

    if np.allclose(rotation, np.eye(3), atol=1e-8):
        return atoms, False

    rotated = atoms.copy()
    rotated.calc = atoms.calc
    rotated.set_cell(np.asarray(atoms.get_cell(), float) @ rotation.T, scale_atoms=True)
    if verbose:
        print(f"[ElasticASE] Reorientation: rotated the structure into the standard setting "
              f"(space group #{number}); cell and atoms unchanged.")
    return rotated, True


class ElasticASE:
    """
    Second-order elastic constants of a crystal from energy-strain calculations.

    Also available as ``ElasticConstants``. The energies come from

    * the attached ASE calculator: ``calculate()`` (in memory) or ``calculate(workdir=...)`` /
      ``run_pending(workdir)`` (one folder per run, energies saved, restarts);
    * run folders of a DFT code: ``write_inputs(workdir)``, then ``run_jobs(workdir)`` or your
      own job scripts, then ``read_results(workdir)``.

    Every route ends in the same :class:`~elastic_ase.results.ElasticResults`, kept as
    ``.results``.
    """

    def __init__(
        self,
        atoms: Atoms,
        symmetry: str = "auto",
        strains=None,
        optimizer=None,
        fmax: float = 1e-3,
        reorient: bool = True,
        symprec: float = 1e-3,
        frame: str = "cartesian",
        verbose: bool = True,
        code: str = None,
        command: str = None,
        output: str = None,
    ):
        """
        Parameters
        ----------
        atoms : ase.Atoms
            Relaxed structure. A calculator is needed for calculate() / run_pending(), not for
            writing or reading run folders.
        symmetry : str
            "auto" (detect with spglib) or a named table: "cubic", "hexagonal",
            "tetragonal_1", "tetragonal_2", "trigonal_1", "trigonal_2", "orthogonal"
            (= "orthorhombic"), "monoclinic", "triclinic".
        strains : list of float, optional
            Strain values (at least 3 different ones); default -1 % to +1 % in steps of 0.5 %.
        optimizer : ASE optimizer class, optional
            Relaxes the atoms inside every strained cell (calculate / run_pending only).
        fmax : float
            Force criterion of the optimizer (eV/Å).
        reorient : bool
            Rotate the structure rigidly into the standard setting (cell and atoms unchanged).
        symprec : float
            spglib tolerance (Å).
        frame : str
            "cartesian" or "standard": frame of ``results["stiffness_matrix"]``.
        verbose : bool
            Print progress and the report (messages also go to the "elastic_ase" logger).
        code : str, optional
            DFT code of the run folders, "vasp" or "qe" (Quantum ESPRESSO). Sets what
            write_inputs() writes, which output file status() / read_results() read and check,
            and the default command of run_jobs().
        command : str, optional
            Shell command run in every run folder by run_jobs(), e.g. "mpirun -np 16 vasp_std".
        output : str, optional
            Output file to read in every run folder; overrides the code's default (e.g.
            "vasprun.xml"), or names the file for another code ASE can read.
        """
        self.verbose = verbose
        self._code = codes.get_code(code)
        self.command = command
        self.output = output
        self.optimizer = optimizer
        self.fmax = fmax
        self.symprec = symprec
        self.frame = {"input": "cartesian", "raw": "cartesian"}.get(frame.lower().strip(), frame.lower().strip())
        if self.frame not in ("cartesian", "standard"):
            raise ValueError(f"frame must be 'cartesian' or 'standard', got {frame!r}")

        self.strains = np.array(sorted(float(s) for s in (DEFAULT_STRAINS if strains is None else strains)))
        if len(np.unique(self.strains)) < 3:
            raise ValueError(
                f"At least 3 different strain values are needed for the quadratic fit, got {self.strains.tolist()}"
            )

        self.atoms = atoms.copy()
        self.atoms.calc = atoms.calc
        self._setup_symmetry(symmetry, reorient)
        self.V0 = self.atoms.get_volume()
        self.results = None
        self.exit_codes = {}

    def _setup_symmetry(self, symmetry, reorient):
        """Space group, rigid rotation into the standard setting, deformation table, and the
        orientation check (all from the structure, before any calculation)."""
        cell = (np.asarray(self.atoms.get_cell(), float), self.atoms.get_scaled_positions(), self.atoms.numbers)
        dataset = spg_call("get_symmetry_dataset", cell, symprec=self.symprec)
        get = (lambda key: getattr(dataset, key)) if hasattr(dataset, "number") else \
            (lambda key: dataset[key])
        self.space_group_number = int(get("number")) if dataset is not None else 0
        self.space_group_symbol = str(get("international")) if dataset is not None else "Unknown"

        self.reoriented = False
        if reorient:
            self.atoms, self.reoriented = standardize_and_reorient(self.atoms, symprec=self.symprec,
                                                                   verbose=self.verbose)
        # rotation from the calculation frame (after any reorientation) to the standard setting
        try:
            self.standard_rotation, _ = auto.standard_setting(self.atoms, symprec=self.symprec)
        except Exception:
            self.standard_rotation = None

        requested = symmetry.lower().strip()
        self.is_auto_oriented = False
        if requested == "auto":
            self.symmetry = auto_detect_symmetry(self.atoms, symprec=self.symprec)
            if reorient:
                self.crystal_system = SYSTEMS[self.symmetry]
                self._say(f"[ElasticASE] Auto-detected crystal symmetry: '{self.symmetry}'")
            else:
                # orientation-aware mode: constants and deformations derived in this orientation
                self.crystal_system, self.basis = auto.from_structure(self.atoms, symprec=self.symprec)
                self.is_auto_oriented = True
                self._say(f"[ElasticASE] Orientation-aware automatic mode: derived {len(self.basis)} "
                          "independent constants in the current Cartesian frame.")
        else:
            self.symmetry = SYSTEM_ALIASES.get(requested, requested)
            if self.symmetry not in SYSTEMS:
                raise ValueError(f"Unknown symmetry: '{symmetry}'. Supported: {sorted(SYSTEMS)} or 'auto'")
            self.crystal_system = SYSTEMS[self.symmetry]

        self.orientation_warnings = [] if self.is_auto_oriented else validate_orientation(
            self.atoms, self.symmetry, verbose=self.verbose, symprec=self.symprec)

    # ----------------------------------------------------------------- helpers
    def _say(self, message):
        """Progress message: to the "elastic_ase" logger, and printed when verbose."""
        logger.info(message)
        if self.verbose:
            print(message)

    def _warn(self, message):
        warnings.warn(message)
        if self.verbose:
            print(f"  [Warning] {message}")

    @property
    def code(self):
        """Name of the DFT code of the run folders ("vasp", "qe") or None."""
        return self._code.name if self._code else None

    def _code_for(self, code=None):
        return codes.get_code(code) or self._code

    def _output_name(self, filename=None, code=None):
        """Output file to read: explicit filename, else this object's output, else the code's."""
        if filename:
            return filename
        adapter = self._code_for(code)
        if self.output and (code is None or adapter is self._code):
            return self.output
        return adapter.output if adapter else None

    def _quietly(self, function, *args, **kwargs):
        verbose, self.verbose = self.verbose, False
        try:
            return function(*args, **kwargs)
        finally:
            self.verbose = verbose

    # ------------------------------------------------------------------- runs
    def runs(self):
        """List of (run name, atoms) for every calculation needed.

        The run name is also the run's folder: ``ground_state`` for the unstrained cell and
        ``<deformation>/strain_<value>`` for each strained cell. A zero strain is the unstrained
        cell, so it is calculated once (as ``ground_state``) and not for every deformation.
        """
        cell = np.asarray(self.atoms.get_cell(), float)
        out = [(GROUND_STATE, self.atoms.copy())]
        for deformation in self.crystal_system.deformations:
            for s in self.strains:
                if s == 0.0:
                    continue
                strained = self.atoms.copy()
                strained.set_cell(cell @ build_strain_matrix(*deformation.at(float(s))).T, scale_atoms=True)
                out.append((strained_key(deformation.folder, float(s)), strained))
        return out

    def _print_header(self):
        self._say("=" * 70)
        self._say(f" ElasticASE: Elastic Constants Calculation ({self.symmetry.upper()})")
        self._say("=" * 70)
        self._say(f"Initial Cell Volume V0 : {self.V0:.4f} Å³ ({self.V0 / len(self.atoms):.4f} Å³/atom)")
        self._say(f"Strain values          : {self.strains.tolist()}")
        self._say(f"Deformation modes      : {len(self.crystal_system.deformations)}")
        self._say(f"Calculations           : {len(self.runs())} (including the unstrained cell)")
        self._say("-" * 70)

    # -------------------------------------------------------------- run state
    def status(self, workdir, filename=None, code=None):
        """
        State of every calculation: {"done": [...], "failed": {run: reason}, "pending": [...]}.

        Read from the output files (``filename``, or the output of ``code`` / of the code given
        to the constructor, e.g. OUTCAR for VASP): each is read and checked (readable, expected
        cell, converged); a missing output counts as pending. Without any of them, from the
        energies saved by calculate(workdir=...) / run_pending().
        """
        filename = self._output_name(filename, code)
        state = rundir.status(workdir, self.runs(), filename, self._code_for(code))
        self._say(f"[ElasticASE] {workdir}: {len(state['done'])} done, {len(state['failed'])} failed, "
                  f"{len(state['pending'])} pending (of {len(self.runs())})")
        for key, reason in state["failed"].items():
            self._say(f"  failed   {key}: {reason}")
        for key in state["pending"]:
            self._say(f"  pending  {key}")
        return state

    def missing(self, workdir, filename=None, code=None):
        """Runs without a usable result yet (failed and pending), in run order.

        Pass them to run_pending(workdir, runs=...) or write_inputs(workdir, runs=...).
        """
        state = self._quietly(self.status, workdir, filename, code)
        failed = set(state["failed"])
        return [key for key, _ in self.runs() if key in failed or key in state["pending"]]

    def is_complete(self, workdir, filename=None, code=None):
        """True when every calculation has a usable result.

        Typical use with a calculator that sometimes fails::

            for attempt in range(3):
                if ec.is_complete(workdir):
                    break
                ec.run_pending(workdir)
            results = ec.calculate(workdir=workdir)
        """
        return not self.missing(workdir, filename, code)

    # ------------------------------------------------- running with the calculator
    def _energy(self, atoms, key, calc):
        """Energy of one run with the attached calculator; raises on failure or no convergence."""
        atoms.calc = calc
        if self.optimizer is not None and key != GROUND_STATE:
            self.optimizer(atoms, logfile=None).run(fmax=self.fmax)
        energy = float(atoms.get_potential_energy())
        if getattr(calc, "converged", None) is False:
            raise RuntimeError("calculator reports converged = False (unconverged SCF)")
        return energy

    def _calculator(self):
        calc = self.atoms.calc
        if calc is None:
            raise ValueError("Attach a calculator first (atoms.calc = ...), or use write_inputs() / read_results().")
        return calc

    def run_pending(self, workdir, runs=None):
        """
        Run the calculations that are not done yet (failed and pending ones), or only ``runs``
        (a list of run names as given by status()), with the attached calculator.

        Each run gets its own folder workdir/<run> (the calculator's ``directory`` is set to
        it, so VASP, Quantum ESPRESSO, ... keep every run's output). Each energy is saved to
        workdir/energies.json as soon as it is computed. A run that raises an error, or whose
        calculator reports ``converged = False`` (ASE's VASP calculator after an unconverged
        SCF), is recorded as failed and the other runs still go ahead.
        Returns the state from the saved energies (see status()).
        """
        calc = self._calculator()
        os.makedirs(workdir, exist_ok=True)
        all_runs = self.runs()
        done = rundir.load_energies(workdir, all_runs)
        if runs is None:
            todo = [(k, a) for k, a in all_runs if k not in done]
        else:
            selected = set(runs)
            unknown = selected - {key for key, _ in all_runs}
            if unknown:
                raise ValueError(f"Unknown runs: {sorted(unknown)}")
            todo = [(k, a) for k, a in all_runs if k in selected]

        has_directory = hasattr(calc, "directory")
        original_directory = getattr(calc, "directory", None)
        try:
            for key, atoms in todo:
                if has_directory:
                    run_dir = os.path.join(workdir, key)
                    os.makedirs(run_dir, exist_ok=True)
                    rundir.set_directory(calc, run_dir)
                try:
                    energy = self._energy(atoms, key, calc)
                except Exception as exc:
                    error = f"{type(exc).__name__}: {exc}"
                    rundir.record(workdir, key, error=error)
                    self._say(f"  [Failed] {key}: {error}")
                    continue
                rundir.record(workdir, key, atoms=atoms, energy=energy)
                self._say(f"  [Done]   {key}: E = {energy:.6f} eV")
        finally:
            if has_directory:
                calc.directory = original_directory
        return rundir.status(workdir, all_runs)

    def calculate(self, workdir=None):
        """
        Run all strained cells with the attached calculator, fit and analyse.

        Without workdir everything stays in memory. With workdir the runs are done by
        run_pending(): one folder per run and every energy saved as soon as it is computed,
        so calling calculate(workdir=...) again (or run_pending) only repeats failed or
        missing runs. If any run failed, a RuntimeError lists them; see status(workdir).
        """
        calc = self._calculator()
        self._print_header()
        n_runs = len(self.runs())

        if workdir is not None:
            state = self.run_pending(workdir)
            if state["failed"] or state["pending"]:
                lines = "\n".join(f"  {k}: {v}" for k, v in state["failed"].items())
                raise RuntimeError(
                    f"{len(state['failed']) + len(state['pending'])} of {n_runs} calculations "
                    f"failed:\n{lines}\nThe successful energies are saved in {workdir}. Fix the problem, "
                    "then call run_pending(workdir) or calculate(workdir=...) again: only these are repeated."
                )
            return self._analyse(rundir.load_energies(workdir, self.runs()))

        energies, failed = {}, {}
        for key, atoms in self.runs():
            try:
                energies[key] = self._energy(atoms, key, calc)
            except Exception as exc:
                failed[key] = f"{type(exc).__name__}: {exc}"
                self._say(f"  [Failed] {key}: {failed[key]}")
        if failed:
            lines = "\n".join(f"  {k}: {v}" for k, v in failed.items())
            raise RuntimeError(
                f"{len(failed)} of {n_runs} calculations failed:\n{lines}\n"
                "Use calculate(workdir=...) to keep the successful energies and rerun only the failed ones."
            )
        return self._analyse(energies)

    # -------------------------------------------------------- run folders (pre / post)
    def write_inputs(self, workdir, runs=None, files=None, code=None, format="vasp", filename=None,
                     calculator_inputs=True):
        """
        Pre-processing: write one folder per calculation (workdir/<run>/), to be run outside
        Python or with run_jobs().

        What is written depends on the code (argument or constructor):

        * ``code="vasp"``: a POSCAR in every folder plus the ``files`` given, e.g.
          ``files={"INCAR": "INCAR", "KPOINTS": "KPOINTS", "POTCAR": "POTCAR"}`` (copied
          unchanged; the POTCAR element order is checked against the POSCAR). With no files
          and an ASE ``Vasp`` calculator attached, the full input is written from it.
        * ``code="qe"``: espresso.pwi written from the attached ASE ``Espresso`` calculator.
        * no code: the attached file-based calculator's input, or else the structure written
          with ``ase.io.write(..., format=format)``; ``files`` are copied into every folder.

        ``runs`` limits this to some runs, e.g. ``runs=ec.missing(workdir)`` to resubmit.
        Returns the list of folders written.
        """
        files = dict(files or {})
        for name, path in files.items():
            if not os.path.isfile(path):
                raise FileNotFoundError(f"{name} file not found: {path}")
        calc = self.atoms.calc
        adapter = self._code_for(code)
        if adapter is not None:
            plan = adapter.plan(self.atoms, calc, files)
        else:
            if filename is None:
                filename = {"vasp": "POSCAR", "espresso-in": "espresso.pwi"}.get(format, f"structure.{format}")
            plan = codes.WritePlan(use_calculator=calculator_inputs and rundir.can_write_input(calc),
                                   format=format, filename=filename)

        selected = None if runs is None else set(runs)
        original_directory = getattr(calc, "directory", None)
        folders = []
        try:
            for key, atoms in self.runs():
                if selected is not None and key not in selected:
                    continue
                folder = os.path.join(workdir, key)
                os.makedirs(folder, exist_ok=True)
                if plan.use_calculator:
                    rundir.set_directory(calc, folder)
                    rundir.write_calculator_input(calc, atoms)
                else:
                    ase_write(os.path.join(folder, plan.filename), atoms, format=plan.format)
                for name, path in files.items():
                    shutil.copyfile(path, os.path.join(folder, name))
                folders.append(folder)
        finally:
            if plan.use_calculator:
                calc.directory = original_directory

        rundir.write_json(workdir, rundir.RUNS_FILE, {
            "symmetry": self.symmetry, "strains": self.strains.tolist(),
            "runs": [key for key, _ in self.runs()]})
        what = "calculator inputs" if plan.use_calculator else f"{plan.filename} files"
        self._say(f"[ElasticASE] Wrote {what} for {len(folders)} calculations to {workdir}")
        for note in plan.notes:
            self._say(f"[ElasticASE] Note: {note}")
        if self.optimizer is not None:
            self._say("[ElasticASE] Note: relax the atoms inside the fixed strained cells in your own "
                      "calculations (e.g. VASP ISIF = 2); the optimizer is not used here.")
        return folders

    def run_jobs(self, workdir, command=None, runs=None, code=None, filename=None):
        """
        Run the DFT code in every run folder that has no usable output yet (or only ``runs``),
        one after the other, then return status().

        ``command`` is run with the shell inside each folder, e.g. ``"mpirun -np 16 vasp_std"``
        or ``"mpirun pw.x -in espresso.pwi > espresso.pwo"``. Default: the constructor's
        command, else the code's default (``vasp_std``, ``pw.x -in espresso.pwi > espresso.pwo``).
        The command's screen output goes to <folder>/elastic_ase_run.log and the exit codes are
        kept in ``.exit_codes``. A run that exits with an error or leaves no usable output shows
        up as failed or pending in the returned status, so simply call run_jobs again after
        fixing the cause::

            ec.write_inputs(workdir, files=...)
            for attempt in range(3):
                if ec.is_complete(workdir):
                    break
                ec.run_jobs(workdir)
            results = ec.read_results(workdir)

        For many runs on a cluster, submit the folders as separate jobs instead and use
        status() / read_results() afterwards.
        """
        adapter = self._code_for(code)
        output = self._output_name(filename, code)
        if output is None:
            raise ValueError("Which output file shows that a run is done? Give code= or filename= "
                             "(here or in the constructor).")
        command = command or self.command or (adapter.command if adapter else None)
        if not command:
            raise ValueError("Give the command to run in every folder (command=...).")

        todo = list(runs) if runs is not None else self.missing(workdir, filename=output, code=code)
        unknown = [key for key in todo if key not in {k for k, _ in self.runs()}]
        if unknown:
            raise ValueError(f"Unknown runs: {unknown}")
        self.exit_codes = {}
        for key in todo:
            folder = os.path.join(workdir, key)
            if not os.path.isdir(folder):
                self.exit_codes[key] = None
                self._say(f"  [Skipped] {key}: no input folder, call write_inputs first")
                continue
            self._say(f"  [Running] {key}: {command}")
            with open(os.path.join(folder, "elastic_ase_run.log"), "w") as log:
                proc = subprocess.run(command, shell=True, cwd=folder, stdout=log, stderr=subprocess.STDOUT)
            self.exit_codes[key] = proc.returncode
            if proc.returncode != 0:
                self._say(f"  [Failed]  {key}: exit code {proc.returncode} (see {folder}/elastic_ase_run.log)")
        return self.status(workdir, filename=output, code=code)

    def read_results(self, workdir, filename=None, code=None):
        """
        Post-processing: read the energy of every calculation from ``workdir/<run>/<filename>``
        (the code's output, e.g. OUTCAR; or any file ase.io.read understands), then fit and
        analyse.

        Each file must hold the strained cell this run expects (catches mixed-up folders), and
        the code's convergence check (VASP, Quantum ESPRESSO) must pass. Missing, unreadable,
        mismatched or unconverged runs are all listed in one RuntimeError; see status() and
        write_inputs(..., runs=...) to resubmit them.
        """
        filename = self._output_name(filename, code) or "OUTCAR"
        self._print_header()
        energies, problems = rundir.read_outputs(workdir, self.runs(), filename, self._code_for(code))
        if problems:
            lines = "\n".join(f"  {k}: {v}" for k, v in problems.items())
            raise RuntimeError(f"{len(problems)} of {len(self.runs())} calculations are not usable:\n{lines}")
        return self._analyse(energies)

    # --------------------------------------------------------- code shortcuts
    def vasp_input(self, workdir, incar=None, kpoints=None, potcar=None, runs=None):
        """VASP run folders: POSCAR plus your INCAR, KPOINTS and POTCAR (copied unchanged).
        Same as write_inputs(workdir, files={...}, code="vasp"). Without files and with an ASE
        ``Vasp`` calculator attached, the full input is written from the calculator."""
        files = {name: path for name, path in (("INCAR", incar), ("KPOINTS", kpoints), ("POTCAR", potcar))
                 if path is not None}
        return self.write_inputs(workdir, runs=runs, files=files, code="vasp")

    def vasp_status(self, workdir, filename=None):
        """status() of the VASP runs in workdir, from their OUTCARs."""
        return self.status(workdir, filename=filename, code="vasp")

    def vasp_output(self, workdir, filename=None):
        """read_results() of the VASP runs in workdir (OUTCAR, or ``filename``)."""
        return self.read_results(workdir, filename=filename, code="vasp")

    def qe_input(self, workdir, runs=None):
        """Quantum ESPRESSO run folders (espresso.pwi) from the attached ASE ``Espresso``
        calculator. Same as write_inputs(workdir, code="qe")."""
        return self.write_inputs(workdir, runs=runs, code="qe")

    def qe_status(self, workdir, filename=None):
        """status() of the Quantum ESPRESSO runs in workdir, from their pw.x outputs."""
        return self.status(workdir, filename=filename, code="qe")

    def qe_output(self, workdir, filename=None):
        """read_results() of the Quantum ESPRESSO runs in workdir (espresso.pwo)."""
        return self.read_results(workdir, filename=filename, code="qe")

    # --------------------------------------------------------------- analysis
    def _analyse(self, energies):
        """Fit, solve and derive everything; keep and return the ElasticResults."""
        results = analysis.analyse(self.crystal_system, self.strains, energies, self.atoms, self.symmetry,
                                   self.standard_rotation, self.space_group_number, self.frame)
        for name, a in results.coefficients.items():
            self._say(f"  Mode {name:<18} -> Coeff = {a:10.2f} GPa  (R² = {results.r2_scores[name]:.6f})")
        message = analysis.residual_stress_message(results.residual_stress)
        if message:
            self._warn(message)
        self.results = results
        self.directional_analysis()  # fills results.directional_properties (needs elastic-vis)
        if self.verbose:
            print(report.format_report(results, self.frame, f"{self.space_group_symbol} (#{self.space_group_number})"))
        return results

    def _print_report(self):
        print(report.format_report(self.results, self.frame, f"{self.space_group_symbol} (#{self.space_group_number})"))

    def _require_results(self):
        if self.results is None:
            raise RuntimeError("No results yet: run calculate(), read_results() or vasp_output() first.")
        return self.results

    # read-only views of the results (kept for compatibility)
    def _result(name, doc):
        return property(lambda self: getattr(self.results, name) if self.results is not None else None, doc=doc)

    stiffness = _result("stiffness_cartesian", "C in the calculation frame (None before an analysis).")
    stiffness_standard = _result("stiffness_standard", "C in the standard setting.")
    constants = _result("independent_constants", "Independent constants of the table.")
    standard_constants = _result("standard_constants", "Independent constants in the standard setting.")
    averages = _result("polycrystalline_averages", "Voigt-Reuss-Hill averages.")
    properties = _result("derived_properties", "Derived properties (density, sound velocities, ...).")
    directional = _result("directional_properties", "Directional extremes (elastic-vis).")
    is_stable = _result("is_stable", "Born stability.")
    stability_criteria = _result("stability_criteria", "Description of the stability criterion.")
    coeffs = _result("coefficients", "Fitted curvature of every deformation (GPa).")
    r2_scores = _result("r2_scores", "R^2 of every fit.")
    fits = _result("fits", "Fit data per deformation.")
    residual_stress = _result("residual_stress", "Residual stress of the reference cell from the fits.")
    energies = _result("energies", "Every energy (eV) by run name.")
    del _result

    @property
    def E0(self):
        """Energy of the unstrained cell (eV)."""
        return self.results.energies[GROUND_STATE] if self.results is not None else None

    @property
    def specs(self):
        """[(folder, deformation)] of the table (kept for compatibility)."""
        return [(d.folder, d.at) for d in self.crystal_system.deformations]

    @property
    def _active_matrix(self):
        r = self._require_results()
        return r.stiffness_standard if self.frame == "standard" else r.stiffness_cartesian

    # ----------------------------------------------------------- plots, files
    def plot_fits(self, filename="energy_strain.png", show=False):
        """
        Plot every energy-strain curve with its quadratic fit (one panel per deformation).
        Smooth parabolas through all points are the first thing to check; scatter, a kink or a
        negative curvature means a problem (unconverged settings, a failed run, too large
        strains). Returns the file name.
        """
        r = self._require_results()
        return plotting.plot_fits(r.fits, r.coefficients, filename, show)

    def save(self, filename="elastic_results.json"):
        """
        Save everything needed to redo the analysis (structure, symmetry, strains and every
        energy) together with the results to a JSON file. ``ElasticConstants.load(filename)``
        recomputes all results from it without running any calculation.
        """
        r = self._require_results()
        from ase.io.jsonio import encode

        data = {
            "format": "elastic_ase results 1",
            "elastic_ase_version": __import__("elastic_ase").__version__,
            "atoms": json.loads(encode(self.atoms.copy())),  # after any reorientation
            "symmetry": self.symmetry,
            "auto_oriented": bool(self.is_auto_oriented),
            "strains": self.strains.tolist(),
            "frame": self.frame,
            "symprec": self.symprec,
            "energies_eV": r.energies,
            "results": _plain({k: v for k, v in r.items() if k not in ("energies", "fits")}),
        }
        with open(filename, "w") as f:
            json.dump(data, f, indent=1)
        self._say(f"[ElasticASE] Saved structure, energies and results to {filename}")
        return filename

    @classmethod
    def load(cls, filename, verbose=False):
        """Recreate from save() and redo the full analysis from the saved energies (no
        calculation is run). The results are in ``.results``."""
        from ase.io.jsonio import decode
        with open(filename) as f:
            data = json.load(f)
        if data.get("format") != "elastic_ase results 1":
            raise ValueError(f"{filename} is not an elastic_ase results file")
        ec = cls(decode(json.dumps(data["atoms"])), symmetry="auto" if data["auto_oriented"] else data["symmetry"],
                 strains=data["strains"], reorient=False, frame=data["frame"], symprec=data["symprec"],
                 verbose=verbose)
        ec._analyse({k: float(v) for k, v in data["energies_eV"].items()})
        return ec

    def export_matrix(self, filepath: str):
        """Write the 6x6 stiffness matrix (frame as chosen) as a plain text file
        (numpy.loadtxt, elastic-vis)."""
        np.savetxt(filepath, self._active_matrix, fmt="%.6f")
        self._say(f"[ElasticASE] Exported 6x6 stiffness matrix to {filepath}")

    @property
    def compliance(self):
        """6x6 Voigt compliance matrix S = C^-1 in 1/GPa (frame as chosen)."""
        return np.linalg.inv(self._active_matrix)

    # ------------------------------------------------------------ elastic-vis
    def to_elastic_vis(self):
        """The stiffness matrix (frame as chosen) as an elastic_vis.Elastic object."""
        from elastic_vis import Elastic
        return Elastic(self._active_matrix)

    def directional_analysis(self):
        """Extremes of Young's modulus, shear modulus, Poisson's ratio and linear
        compressibility over all directions (needs elastic-vis; None without it). Also stored
        in ``results.directional_properties``."""
        r = self._require_results()
        r.directional_properties = self._directional_extremes()
        return r.directional_properties

    def _directional_extremes(self):
        try:
            material = self.to_elastic_vis()
            from elastic_vis.core import minimize, maximize
        except Exception:
            return None
        try:
            _, e_min = minimize(lambda x: material.Young(x), 2)
            _, e_max = maximize(lambda x: material.Young(x), 2)
            _, g_min = minimize(lambda x: material.shear(x), 3)
            _, g_max = maximize(lambda x: material.shear(x), 3)
            _, nu_min = minimize(lambda x: material.Poisson(x), 3)
            _, nu_max = maximize(lambda x: material.Poisson(x), 3)
            _, lc_min = minimize(lambda x: material.LC(x), 2)
            _, lc_max = maximize(lambda x: material.LC(x), 2)
        except Exception as exc:
            logger.info("directional analysis skipped: %s", exc)
            return None
        return {
            "E_min": float(e_min), "E_max": float(e_max),
            "E_anisotropy": float(e_max / e_min) if e_min > 1e-6 else float("nan"),
            "G_min": float(g_min), "G_max": float(g_max),
            "G_anisotropy": float(g_max / g_min) if g_min > 1e-6 else float("nan"),
            "nu_min": float(nu_min), "nu_max": float(nu_max), "is_auxetic": bool(nu_min < 0.0),
            "LC_min": float(lc_min), "LC_max": float(lc_max), "is_negative_LC": bool(lc_min < 0.0),
        }

    def Young(self, theta: float, phi: float) -> float:
        """Young's modulus along the direction (theta, phi) in GPa (elastic-vis)."""
        return self.to_elastic_vis().Young([theta, phi])

    def linear_compressibility(self, theta: float, phi: float) -> float:
        """Linear compressibility along (theta, phi) in TPa^-1 (elastic-vis)."""
        return self.to_elastic_vis().LC([theta, phi])

    def shear(self, theta: float, phi: float) -> tuple:
        """(min, max) shear modulus for the direction (theta, phi) in GPa (elastic-vis)."""
        return self.to_elastic_vis().shear2D([theta, phi])

    def poisson(self, theta: float, phi: float) -> tuple:
        """(min, max) Poisson's ratio for the direction (theta, phi) (elastic-vis)."""
        return self.to_elastic_vis().Poisson2D([theta, phi])

    def visualize(self, properties=("Young", "Shear"), dimensions=("2D", "3D"), prefix="elastic",
                  style_config=None, show=False, output_format="png"):
        """
        Directional plots with elastic-vis.

        properties: "Young", "Shear", "Poisson", "LC" (or "all"); dimensions: "2D", "3D" (or
        "both"); prefix: file name prefix; style_config: dict or path to a style JSON file;
        output_format: "png", "pdf" or "svg". Returns the list of files written.
        """
        import elastic_vis
        if isinstance(style_config, str):
            with open(style_config) as f:
                style_config = json.load(f)
        if isinstance(properties, str):
            properties = ["Young", "Shear", "Poisson", "LC"] if properties.lower() == "all" else [properties]
        if isinstance(dimensions, str):
            dimensions = ["2D", "3D"] if dimensions.lower() in ("both", "all") else [dimensions]
        wanted = {p.lower() for p in properties}
        if "linear_compressibility" in wanted:
            wanted.add("lc")
        dims = {d.upper() for d in dimensions}

        material = elastic_vis.Elastic(self._active_matrix)
        files = []
        for key, plot_name, file_stem in _VIS_PLOTS:
            prop, dim = key
            if prop in wanted and dim in dims:
                getattr(elastic_vis, plot_name)(material, prefix=prefix, style_config=style_config,
                                                show=show, output_format=output_format)
                files.append(f"{prefix}_{file_stem}.{output_format}")
        self._say(f"[ElasticASE] Visualization saved: {', '.join(files)}")
        return files

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

    def __repr__(self):
        status = "calculated" if self.results is not None else "uncalculated"
        return f"<ElasticASE(symmetry='{self.symmetry}', status='{status}')>"

    def __str__(self):
        if self.results is None:
            return f"ElasticASE(symmetry='{self.symmetry}', status='uncalculated')"
        stable = "Stable" if self.results.is_stable else "Unstable"
        return f"ElasticASE(symmetry='{self.symmetry}', C11={self._active_matrix[0, 0]:.2f} GPa, {stable})"


# (property, dimension) -> elastic-vis plot function and the file name it writes
_VIS_PLOTS = [
    (("young", "2D"), "plot_youngs_2d", "Youngs_Modulus_Polar_2D"),
    (("young", "3D"), "plot_youngs_3d", "Youngs_Modulus_3D"),
    (("shear", "2D"), "plot_shear_2d", "Shear_Modulus_Polar_2D"),
    (("shear", "3D"), "plot_shear_3d", "Shear_Modulus_3D"),
    (("poisson", "2D"), "plot_poisson_2d", "Poisson_Ratio_Polar_2D"),
    (("poisson", "3D"), "plot_poisson_3d", "Poisson_Ratio_3D"),
    (("lc", "2D"), "plot_lc_2d", "Linear_Compressibility_Polar_2D"),
    (("lc", "3D"), "plot_lc_3d", "Linear_Compressibility_3D"),
]


def _plain(value):
    """JSON-friendly copy (numpy arrays and scalars to lists and Python numbers)."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (np.ndarray, list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    return value
