"""Run folders: names, saved energies, failure records, and reading outputs back.

Layout of a work directory::

    workdir/
      ground_state/                 unstrained cell
      <deformation>/strain_<value>/ one per deformation and strain
      energies.json                 energies of live runs
      failed.json                   reasons of failed live runs
      elastic_ase_runs.json         symmetry, strains, run names
"""
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
    """Write atomically, so a killed job never leaves a half-written file."""
    path = os.path.join(workdir, name)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=1)
    os.replace(tmp, path)


def load_energies(workdir, runs):
    """Saved energies in workdir/energies.json that belong to these runs ([(key, atoms)]):
    an entry is used only if its cell and number of atoms match the run."""
    expected = dict(runs)
    energies = {}
    for key, entry in read_json(workdir, ENERGY_FILE).items():
        if key in expected and same_cell(entry["cell"], expected[key].get_cell()) \
                and entry.get("natoms") == len(expected[key]):
            energies[key] = float(entry["energy"])
    return energies


def record(workdir, key, atoms=None, energy=None, error=None):
    """Save a finished run (energy) or a failed one (error)."""
    saved = read_json(workdir, ENERGY_FILE)
    failed = read_json(workdir, FAILED_FILE)
    if error is None:
        saved[key] = {"energy": float(energy), "cell": np.asarray(atoms.get_cell()).tolist(),
                      "natoms": len(atoms)}
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
        reason = codes.not_converged(path, code)
        if reason:
            problems[key] = reason
            continue
        energies[key] = float(energy)
    return energies, problems


def status(workdir, runs, filename=None, code=None):
    """{"done": [...], "failed": {run: reason}, "pending": [...]} for these runs.

    From the output files when ``filename`` is given (missing output = pending, unusable =
    failed); otherwise from the energies and failures saved by live runs.
    """
    keys = [key for key, _ in runs]
    if filename is None:
        done = load_energies(workdir, runs)
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
    """Point a file-based calculator at another run folder, keeping its type (str or Path)."""
    calc.directory = Path(path) if isinstance(calc.directory, Path) else str(path)


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
