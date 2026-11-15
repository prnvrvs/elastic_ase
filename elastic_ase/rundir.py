"""Run folders: names, saved energies, failure records, and reading outputs back.

Layout of a work directory::

    workdir/
      ground_state/                 unstrained cell
      <deformation>/strain_<value>/ one per deformation and strain
      energies.json                 energies of live runs
      failed.json                   reasons of failed live runs
      elastic_ase_runs.json         symmetry, strains, run names
"""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from ase.io import read as ase_read

from elastic_ase import codes

GROUND_STATE = "ground_state"
ENERGY_FILE = "energies.json"
FAILED_FILE = "failed.json"
RUNS_FILE = "elastic_ase_runs.json"


def strained_key(mode, strain):
    """Run name of one strained calculation, e.g. 'C11_C12_I/strain_+0.005000'."""
    return f"{mode}/strain_{strain:+.6f}"


def same_cell(cell_a, cell_b, tol=1e-4):
    """True if two cells agree to a relative tolerance (used to match saved runs to runs)."""
    a, b = np.asarray(cell_a, float), np.asarray(cell_b, float)
    return a.shape == b.shape and np.abs(a - b).max() <= tol * max(1.0, np.abs(b).max())


def read_json(workdir, name):
    path = os.path.join(workdir, name)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def write_json(workdir, name, data):
    """Write atomically, so a killed job never leaves a half-written file. Under MPI only
    rank 0 writes; the other ranks wait until the file is there."""
    from ase.parallel import world
    path = os.path.join(workdir, name)
    if world.rank == 0:
        tmp = f"{path}.{os.getpid()}.tmp"
        with open(tmp, "w") as f:
            json.dump(data, f, indent=1)
        os.replace(tmp, path)
    world.barrier()


def fingerprint(atoms):
    """Short hash of the atoms (elements and scaled positions) of a run before any relaxation,
    so saved energies are not reused for a different structure in the same cell."""
    data = np.round(np.asarray(atoms.get_scaled_positions(wrap=True), float), 6)
    data = np.where(data >= 1.0, 0.0, data) + 0.0  # 1.0 and -0.0 after rounding
    digest = hashlib.sha1(np.asarray(atoms.numbers, np.int64).tobytes() + data.tobytes())
    return digest.hexdigest()[:16]


NO_RELAXATION = {"relaxation": False}


def load_energies(workdir, runs, settings=None):
    """Saved energies in workdir/energies.json that belong to these runs ([(key, atoms)]).

    An entry is used only if its cell, atoms (elements and positions before relaxation) and
    relaxation settings match the run; entries saved by older versions carry no atoms or
    settings and count as unrelaxed runs."""
    settings = NO_RELAXATION if settings is None else settings
    expected = dict(runs)
    energies = {}
    for key, entry in read_json(workdir, ENERGY_FILE).items():
        if key not in expected:
            continue
        atoms = expected[key]
        if not same_cell(entry["cell"], atoms.get_cell()) or entry.get("natoms") != len(atoms):
            continue
        if "atoms" in entry and entry["atoms"] != fingerprint(atoms):
            continue
        if entry.get("settings", NO_RELAXATION) != settings:
            continue
        energies[key] = float(entry["energy"])
    return energies


def record(workdir, key, atoms=None, energy=None, error=None, settings=None, structure=None):
    """Save a finished run (energy) or a failed one (error). ``structure`` is the fingerprint
    of the run's atoms before relaxation, ``settings`` the relaxation settings."""
    saved = read_json(workdir, ENERGY_FILE)
    failed = read_json(workdir, FAILED_FILE)
    if error is None:
        saved[key] = {"energy": float(energy), "cell": np.asarray(atoms.get_cell()).tolist(),
                      "natoms": len(atoms),
                      "atoms": structure if structure is not None else fingerprint(atoms),
                      "settings": NO_RELAXATION if settings is None else settings}
        failed.pop(key, None)
        write_json(workdir, ENERGY_FILE, saved)
    else:
        failed[key] = error
    write_json(workdir, FAILED_FILE, failed)


def read_outputs(workdir, runs, filename, code=None):
    """(energies, problems) from workdir/<run>/<filename> for every run ([(key, atoms)]).

    A run is usable if its file exists, ase.io.read gives an energy, the cell and number of atoms
    match the run (catches mixed-up folders) and the code's convergence check does not fail.
    """
    energies, problems = {}, {}
    for key, expected in runs:
        path = os.path.join(workdir, key, filename)
        if not os.path.exists(path):
            problems[key] = "missing " + filename
            continue
        # The convergence check reads the text itself, so it also names an unconverged run whose
        # output ASE cannot parse (pw.x stopped after "convergence NOT achieved").
        reason = codes.not_converged(path, code)
        if reason:
            problems[key] = reason
            continue
        try:
            result = ase_read(path, index=-1)
            energy = result.get_potential_energy()
        except Exception as exc:
            problems[key] = f"cannot read {filename} ({type(exc).__name__}: {exc})"
            continue
        if len(result) != len(expected) or not same_cell(result.get_cell(), expected.get_cell()):
            problems[key] = (f"{filename} holds a different cell than this run expects "
                             "(another structure, strain set or orientation?)")
            continue
        energies[key] = float(energy)
    return energies, problems


def status(workdir, runs, filename=None, code=None, settings=None):
    """{"done": [...], "failed": {run: reason}, "pending": [...]} for these runs.

    From the output files when ``filename`` is given (missing output = pending, unusable =
    failed); otherwise from the energies and failures saved by live runs.
    """
    keys = [key for key, _ in runs]
    if filename is None:
        done = load_energies(workdir, runs, settings)
        saved_failures = read_json(workdir, FAILED_FILE)
        failed = {k: saved_failures[k] for k in keys if k not in done and k in saved_failures}
    else:
        done, problems = read_outputs(workdir, runs, filename, code)
        failed = {k: v for k, v in problems.items() if not v.startswith("missing ")}
    return {
        "done": [k for k in keys if k in done],
        "failed": failed,
        "pending": [k for k in keys if k not in done and k not in failed],
    }


def set_directory(calc, path):
    """Point a file-based calculator at another run folder, keeping its type (str or Path).

    The calculation then runs inside that folder, so a relative pseudopotential folder
    (QE ``pseudo_dir``) is made absolute first: relative to where you started,
    not to the run folder."""
    _absolute_pseudo_dir(calc)
    calc.directory = Path(path) if isinstance(calc.directory, Path) else str(path)


def _absolute_pseudo_dir(calc):
    profile = getattr(calc, "profile", None)
    pseudo = getattr(profile, "pseudo_dir", None)
    if pseudo and not os.path.isabs(str(pseudo)):
        profile.pseudo_dir = os.path.abspath(str(pseudo))
    try:
        control = calc.parameters["input_data"]["control"]
    except (AttributeError, KeyError, TypeError):
        return
    if control.get("pseudo_dir") and not os.path.isabs(str(control["pseudo_dir"])):
        control["pseudo_dir"] = os.path.abspath(str(control["pseudo_dir"]))


def write_calculator_input(calc, atoms):
    """Write the calculator's input files into calc.directory without running it."""
    if hasattr(calc, "write_input"):            # Vasp and other FileIOCalculators
        calc.write_input(atoms)
    else:                                       # GenericFileIOCalculator (Espresso, Abinit, ...)
        Path(calc.directory).mkdir(parents=True, exist_ok=True)
        calc.write_inputfiles(atoms, ["energy"])


def can_write_input(calc):
    """True for a file-based ASE calculator that can write its input without running."""
    return calc is not None and hasattr(calc, "directory") and \
        (hasattr(calc, "write_input") or hasattr(calc, "write_inputfiles"))
