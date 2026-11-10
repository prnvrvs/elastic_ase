# ElasticASE

**ElasticASE** is a lightweight, high-performance Python package for calculating second-order elastic stiffness constants ($C_{ij}$) using the **Atomic Simulation Environment (ASE)**.

It is the direct ASE-native counterpart and modern successor to `elastic_vasp`, operating **100% in-memory** with any ASE calculator (including foundation machine learning potentials such as **GRACE**, **MACE**, **CHGNet**, analytical potentials like **EMT**, or DFT calculators like **GPAW**, **VASP**, and **Quantum ESPRESSO**).

---

## Key Features

- **Works with ANY ASE Calculator:** Simply attach `atoms.calc = calc`.
- **Zero VASP Dependency:** Runs entirely in Python/RAM without dumping dozens of temporary directories or POSCAR/INCAR files to disk.
- **CLI Included:** Run directly from your command line: `elastic-ase POSCAR -c grace -m GRACE-2L-OAM`.
- **Automatic Symmetry Detection:** Automatically determines space group and crystal symmetry using `spglib`.
- **Supports All 7 Crystal Systems:**
  - `cubic` (3 constants: $C_{11}, C_{12}, C_{44}$)
  - `hexagonal` (5 constants: $C_{11}, C_{12}, C_{13}, C_{33}, C_{44}$)
  - `orthogonal` (orthorhombic / tetragonal, up to 9 constants)
  - `trigonal_1` & `trigonal_2` (6 / 7 constants)
  - `monoclinic` (13 constants)
  - `triclinic` (21 constants)
- **Born Mechanical Stability Validation:** Automatically checks theoretical stability criteria.
- **Polycrystalline Moduli:** Computes Voigt, Reuss, and Hill Bulk ($K$), Shear ($G$), Young's ($E$) moduli, Poisson's ratio ($\nu$), Pugh ductility ratio ($K/G$), Cauchy pressure, and Zener anisotropy.
- **Internal Atomic Relaxation:** Supports optional ASE optimizers (e.g., `BFGS`) for non-affine multi-atom unit cells.

---

## Installation

```bash
git clone https://github.com/<your-username>/elastic_ase.git
cd elastic_ase
pip install .
```

For development (editable mode):
```bash
pip install -e .
```

---

## Command-Line Interface (CLI)

Calculate elastic constants directly from any structure file (CIF, POSCAR, XYZ, etc.):

```bash
# Quick test using ASE's built-in EMT calculator:
elastic-ase POSCAR -c emt

# Using GRACE machine learning foundation model:
elastic-ase POSCAR -c grace -m GRACE-2L-OAM

# Save output to JSON or CSV:
elastic-ase structure.cif -c grace -m GRACE-2L-OAM -o results.json
elastic-ase structure.cif -c grace -m GRACE-2L-OAM -o results.csv
```

### CLI Options
- `-c, --calculator`: `emt` or `grace` (default: `emt`).
- `-m, --model`: GRACE model name (default: `GRACE-2L-OAM`).
- `-s, --symmetry`: Crystal symmetry family (`auto`, `cubic`, `hexagonal`, etc., default: `auto`).
- `--strains`: Comma-separated strain list (default: `-0.01,-0.005,0.0,0.005,0.01`).
- `-o, --output`: File path to export results (`.json` or `.csv`).
- `-q, --quiet`: Suppress verbose terminal logs.

---

## Python API Quickstart

### 1. Minimal Example (Instant with EMT)
```python
from ase.build import bulk
from ase.calculators.emt import EMT
from elastic_ase import ElasticASE

# Setup structure and calculator
al = bulk("Al", "fcc", a=4.05, cubic=True)
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

## Running Unit Tests

Run the test suite with standard `unittest`:
```bash
python -m unittest discover -s tests
```
Or with `pytest`:
```bash
pytest tests/
```
