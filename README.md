# elastic_ase

Second-order elastic constants ($C_{ij}$) from energy–strain calculations with ASE
calculators. elastic_ase is the ASE counterpart of
[elastic_vasp](https://github.com/prnvrvs/elastic_vasp): the same deformation tables and
analysis, with the energies from an ASE calculator (run in memory) or from DFT run folders
(VASP, Quantum ESPRESSO) that elastic_ase writes and reads back.

**Documentation:** https://prnvrvs.github.io/elastic_ase/

## Features

- Symmetry detection with spglib and deformation tables for all crystal systems: cubic,
  hexagonal, tetragonal (4/mmm and 4/m), trigonal (with $C_{14}$ and $C_{15}$),
  orthorhombic, monoclinic, triclinic; an orientation-aware automatic mode.
- The structure is rotated rigidly into the standard setting (cell, atoms and per-atom data
  unchanged); without rotation, an exact check warns when a table does not fit the
  orientation.
- Energies from an ASE calculator (`calculate()`), or run folders for a DFT code with
  done / failed / pending tracking and restarts (`code="vasp"` or `"qe"`, `write_inputs`,
  `run_jobs`, `status`, `read_results`).
- Results: stiffness matrix, Born stability, Voigt–Reuss–Hill averages, sound velocities,
  Debye temperature, residual pressure of the reference cell, energy–strain fit plots,
  save/load; directional properties and plots with
  [elastic-vis](https://github.com/prabhat7iitm/elastic-vis).
- Command line tool (`elastic-ase`) for EMT, GRACE or any calculator given as `-c module:function`.

## What is tested

- All tables recover a known stiffness matrix exactly from a synthetic quadratic energy, and
  they are identical to elastic_vasp's.
- With EMT: relaxed alloy structures of every crystal family give the same constants from
  the named table, the full triclinic fit and the automatic mode on randomly rotated copies.
- VASP: reading 19 real VASP 6 OUTCARs of fcc Al reproduces elastic_vasp's result
  ($C_{11}$ = 106.6, $C_{12}$ = 63.4, $C_{44}$ = 28.9 GPa).
- Quantum ESPRESSO: input writing and output reading/convergence checks are tested on files,
  **not yet with a real pw.x run through elastic_ase**.
- **Not tested here:** GRACE, MACE, CHGNet or other ML potentials (they should work like any
  ASE calculator), ASE's `Vasp` calculator writing full inputs, `run_jobs` with a real DFT
  code (tested with a stand-in program), other DFT codes, magnetic systems.

Please report problems on the [issue tracker](https://github.com/prnvrvs/elastic_ase/issues).

---

## Installation

```bash
pip install elastic-ase                                          # from PyPI (once the first release is published)
pip install "git+https://github.com/prnvrvs/elastic_ase.git"     # current version from GitHub
```

Optional, for the directional plots and extremes:
[elastic-vis](https://github.com/prabhat7iitm/elastic-vis) is not on PyPI, so install it
separately:

```bash
pip install "git+https://github.com/prabhat7iitm/elastic-vis.git"
```

For development (editable mode):
```bash
git clone https://github.com/prnvrvs/elastic_ase.git
cd elastic_ase
pip install -e ".[test]"
```

Requires Python 3.10 or newer.

---

## Command-Line Interface (CLI)

Calculate elastic constants directly from any structure file (CIF, POSCAR, XYZ, etc.):

```bash
# Quick test using ASE's built-in EMT calculator:
elastic-ase POSCAR -c emt

# Using GRACE machine learning foundation model:
elastic-ase POSCAR -c grace -m GRACE-2L-OAM

# Symmetry check without running calculations:
elastic-ase POSCAR --check

# Save output to JSON, CSV, or plain 6x6 text matrix:
elastic-ase structure.cif -c grace -m GRACE-2L-OAM -o results.json
elastic-ase structure.cif -c grace -m GRACE-2L-OAM -o results.csv
elastic-ase structure.cif -c grace -m GRACE-2L-OAM -o C.txt

# Calculate and generate 2D & 3D directional plots (via elastic-vis):
elastic-ase structure.cif -c grace -m GRACE-2L-OAM --plot --Young --Shear --2D --3D --format png
elastic-ase structure.cif -c grace -m GRACE-2L-OAM --plot --Poisson --LC --format pdf
```

### CLI Options
- `--check`: Instant symmetry and orientation check without calculations.
- `-c, --calculator`: `emt` (default), `grace`, or any ASE calculator as `module:function`, e.g. `ase.calculators.emt:EMT` (tested with EMT only).
- `-m, --model`: GRACE model name (default: `GRACE-2L-OAM`).
- `-s, --symmetry`: Crystal symmetry family (`auto`, `cubic`, `hexagonal`, `tetragonal_1`, `tetragonal_2`, `orthogonal`, `trigonal_1`, `trigonal_2`, `monoclinic`, `triclinic`).
- `--strains`: Comma-separated strain list (default: `-0.01,-0.005,0.0,0.005,0.01`).
- `-o, --output`: File path to export results (`.json`, `.csv`, or `.txt` for 6x6 matrix).
- `-q, --quiet`: Suppress verbose terminal logs.
- `--plot`: Enable anisotropic property visualization.
- `--Young`, `--Shear`, `--Poisson`, `--LC`: Property selection.
- `--2D`, `--3D`: Dimension selection.
- `--format`: Image format (`png`, `pdf`, `svg`).
- `--style`: Path to custom JSON styling configuration.

---

## Python API Quickstart

### 1. Minimal Example (Instant with EMT)
```python
from ase.build import bulk
from ase.calculators.emt import EMT
from elastic_ase import ElasticASE

# Setup structure and calculator
# Use the relaxed lattice constant of *your* calculator: a pre-stressed cell biases
# the constants. 3.9943 A is the EMT equilibrium of Al (a = 4.05 A would be under
# 1.6 GPa of tension with EMT).
al = bulk("Al", "fcc", a=3.9943, cubic=True)
al.calc = EMT()

# Run ElasticASE (symmetry auto-detected)
el = ElasticASE(al)
res = el.calculate()

print("C11:", res["stiffness_matrix"][0, 0], "GPa")
print("C12:", res["stiffness_matrix"][0, 1], "GPa")
print("C44:", res["stiffness_matrix"][3, 3], "GPa")
print("Bulk Modulus K_VRH:", res["polycrystalline_averages"]["K_VRH"], "GPa")
print("Stable?", res["is_stable"])
```

### 2. With GRACE Foundation Models

Not tested here (GRACE is not installed in our test setup); it should work like any ASE
calculator. Relax the structure with GRACE first.

```python
from ase.build import bulk
from tensorpotential.calculator import grace_fm
from elastic_ase import ElasticASE

calc = grace_fm("GRACE-2L-OAM")
co = bulk("Co", "hcp", a=2.498, c=4.053)
co.calc = calc

el = ElasticASE(co, symmetry="hexagonal")
res = el.calculate()
```

---

## DFT codes: VASP, Quantum ESPRESSO, ... (Python API)

The calculations are independent single-point runs (13 for a cubic crystal with 5 strains;
the unstrained cell is run once). With a DFT code there are two ways to run them.

### A. Let ASE run them, with run folders and restart

The run folders, saved energies and restarts are tested with a stand-in file-based
calculator; running ASE's `Vasp` calculator itself through elastic_ase has not been tested yet.

```python
from ase.calculators.vasp import Vasp

atoms.calc = Vasp(xc="PBE", encut=520, kpts=(12, 12, 12), ismear=1, sigma=0.1,
                  ediff=1e-7, isif=2, ibrion=-1)  # your own converged settings
el = ElasticASE(atoms)
results = el.calculate(workdir="elastic_runs")
```

Each run gets its own folder (`elastic_runs/C11_C12_I/strain_+0.005000/`, ...), so every
OUTCAR is kept, and each energy is saved to `elastic_runs/energies.json` as soon as it is
computed. A run that crashes, or whose SCF did not converge (ASE's `converged = False`), is
recorded as failed; the other runs still go ahead and the failed ones are listed at the end.

To check and repeat only what is missing:

```python
el.status("elastic_runs")        # done / failed (with reason) / pending
el.is_complete("elastic_runs")   # True when every run has a usable energy
el.missing("elastic_runs")       # the failed and pending runs
el.run_pending("elastic_runs")   # run only those (or runs=[...] for a chosen few)
results = el.calculate(workdir="elastic_runs")  # fits once complete; runs nothing already done
```

### B. Run folders for the DFT code itself (e.g. on a cluster)

Tell ElasticConstants once which code you use (and, if it should run the jobs, the command);
every later step then knows what to write, which output to read and check, and what to run:

```python
from ase.io import read
from elastic_ase import ElasticConstants

ec = ElasticConstants(read("POSCAR"), code="vasp", command="mpirun -np 16 vasp_std")
ec.write_inputs("runs", files={"INCAR": "INCAR", "KPOINTS": "KPOINTS", "POTCAR": "POTCAR"})
# -> runs/<deformation>/strain_<value>/{POSCAR, INCAR, KPOINTS, POTCAR}

ec.run_jobs("runs")        # runs the command in every folder without a usable OUTCAR
ec.status("runs")          # done / failed (unconverged, unreadable, wrong cell) / pending
ec.is_complete("runs")     # True when every run has a usable OUTCAR
ec.missing("runs")         # the failed and pending runs
ec.run_jobs("runs")        # call again after fixing a problem: only missing runs are repeated
results = ec.read_results("runs")
```

`run_jobs` runs the folders one after the other (each command's screen output goes to
`elastic_ase_run.log` in its folder). On a cluster, submit the folders as separate jobs
instead, then use the same `status`, `missing` (for `write_inputs(..., runs=...)` to
resubmit) and `read_results`.

| `code=` | writes | reads and checks | default `command` |
|---|---|---|---|
| `"vasp"` | POSCAR + your `files` (or full input from an ASE `Vasp` calculator) | `OUTCAR` (electronic convergence) | `vasp_std` |
| `"qe"` | `espresso.pwi` from an ASE `Espresso` calculator | `espresso.pwo` (SCF convergence) | `pw.x -in espresso.pwi > espresso.pwo` |
| none | structure files (`format=`) + your `files` | `output=` (any file ASE reads) | `command=` |

`output="vasprun.xml"` reads another file than the code's default. The shortcuts
`vasp_input` / `vasp_status` / `vasp_output` and `qe_input` / `qe_status` / `qe_output` do the
same with the code fixed:

```python
ec = ElasticConstants(read("POSCAR"))
ec.vasp_input("runs", incar="INCAR", kpoints="KPOINTS", potcar="POTCAR")
results = ec.vasp_output("runs")
```

Your INCAR, KPOINTS and POTCAR are copied unchanged into every folder (use `ISIF = 2`,
`IBRION = -1` or a fixed-cell ionic relaxation); the POTCAR element order is checked against
the POSCAR. With an ASE `Vasp(...)` calculator attached and no files given, `vasp_input`
writes the full input from the calculator instead.

For Quantum ESPRESSO the input is written from an ASE `Espresso` calculator (its
pseudopotentials, k-points and `input_data` are used for every run):

```python
from ase.calculators.espresso import Espresso, EspressoProfile

atoms.calc = Espresso(profile=EspressoProfile(command="pw.x", pseudo_dir="pseudo"),
                      pseudopotentials={"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
                      input_data={"system": {"ecutwfc": 60, "occupations": "smearing",
                                             "smearing": "mv", "degauss": 0.02}},
                      kpts=(12, 12, 12))
ec = ElasticConstants(atoms, code="qe", command="mpirun -np 16 pw.x -in espresso.pwi > espresso.pwo")
ec.write_inputs("runs")             # runs/<run>/espresso.pwi
ec.run_jobs("runs")                 # or submit the folders yourself
ec.status("runs")                   # done / failed (SCF not converged, unreadable) / pending
results = ec.read_results("runs")   # read, check, fit, analyse
```

Any other code ASE can read: `ElasticConstants(atoms, output="<output file>",
command="<command>")` with `write_inputs(workdir, format=..., files=...)`.

`read_results` checks that every output is readable, holds the strained cell its folder
expects (catches mixed-up folders) and, for VASP, that the electronic loop converged;
everything that is not usable is listed in one error. Use the same structure, symmetry and
strains for writing and reading.

Do your own convergence tests (cutoff, k-points, smearing): converge the elastic constants
themselves, not only the total energy, and use the same settings for every run. Relax the
reference structure first; with `optimizer=` (or `ISIF = 2` in VASP) the atoms are relaxed
inside each fixed strained cell.

---

## Results, checks and derived properties

`calculate()`, `vasp_output()`, `qe_output()` and `read_results()` all return the same
dictionary (also kept as `ec.results`):

| Key | Content |
|---|---|
| `stiffness_matrix` | 6×6 C in GPa (Cartesian frame, or standard setting with `frame="standard"`) |
| `independent_constants`, `standard_constants` | e.g. `C11`, `C12`, `C44` |
| `is_stable`, `stability_criteria` | Born stability (all eigenvalues of C positive) |
| `polycrystalline_averages` | Voigt/Reuss/Hill K, G, E, ν, Pugh ratio, Zener and universal anisotropy |
| `derived_properties` | density, longitudinal/transverse/mean sound velocity, Debye temperature (Anderson), Kleinman parameter (cubic), empirical Vickers hardness (Chen, Tian; rough) |
| `residual_stress` | residual pressure of the reference cell, from the linear terms of the fits |
| `directional_properties` | min/max of E, G, ν, linear compressibility (with elastic-vis) |
| `coefficients`, `r2_scores`, `energies` | fitted curvatures, fit quality, every energy |

Checks and tools:

- **Residual stress:** a reference cell that is not relaxed biases the constants by about the
  same amount. ElasticASE estimates the residual pressure from the linear terms of the
  energy–strain fits (no stress calculation needed, so it also works on OUTCARs) and warns
  above 0.5 GPa. The linear term is taken from a cubic fit, so anharmonicity at large strains
  does not show up as a fake pressure.
- **Look at the curves:** `ec.plot_fits("energy_strain.png")` plots every energy–strain curve
  with its fit. Smooth parabolas are the first check; scatter or kinks point to unconverged
  settings or a failed run.
- **Save and reload:** `ec.save("al.json")` stores the structure, settings and every energy
  with the results; `ElasticASE.load("al.json")` redoes the whole analysis without running
  anything.
- `ElasticConstants` is the same class under a more descriptive name:
  `from elastic_ase import ElasticConstants; ec = ElasticConstants(atoms)`.

---

## Running Unit Tests

Run the test suite with standard `unittest`:
```bash
python -m unittest discover -s tests
```
Or with `pytest`:
```bash
pytest tests/
```

---

## Crystallographic Conventions & Common Pitfalls

Points where the orientation or setting of the crystal matters:

1. **Trigonal Systems & "31m"-Type Groups ($P\bar{3}1m, P\bar{3}1c, P31m, P31c, P312, P3_1 12, P3_2 12$):**
   - In standard settings with $a \parallel x$, the 2-fold axes in 31m-type groups are perpendicular to $a$, causing the stiffness matrix to have $C_{15}$ rather than $C_{14}$.
   - **Solution in ElasticASE:** Automatic symmetry detection routes all trigonal crystals to `trigonal_2`. It fits both $C_{14}$ and $C_{15}$, so it does not depend on the in-plane orientation about $z$ (one extra deformation).

2. **Orientation of the crystal (rotated cells, rhombohedral cells, monoclinic unique axis):**
   - The named tables assume the crystal axes along fixed Cartesian directions: cubic and orthorhombic axes along $x, y, z$; the 3-, 4- or 6-fold axis along $z$; the monoclinic unique axis $b$ along $y$.
   - **Default (`reorient=True`):** the structure is rigidly rotated into this standard setting. Only the orientation changes: your cell (primitive, conventional or supercell), the atoms and per-atom data (magnetic moments, charges, tags, constraints) are kept, so calculator settings such as the k-point mesh stay valid. A primitive rhombohedral cell stays primitive, with its 3-fold axis rotated onto $z$.
   - **`reorient=False` with a named symmetry:** the crystal is used as given. ElasticASE checks exactly whether the chosen table can measure the crystal in this orientation and warns if not (the constants would be wrong).
   - **`reorient=False` with `symmetry="auto"`:** the constants and deformations are derived from the crystal in its own orientation, so any orientation works.

