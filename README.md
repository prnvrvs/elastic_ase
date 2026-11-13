# ElasticASE

**ElasticASE** is a lightweight, high-performance Python package for calculating second-order elastic stiffness constants ($C_{ij}$) using the **Atomic Simulation Environment (ASE)**.

It is the direct ASE-native counterpart and modern successor to `elastic_vasp`, operating **100% in-memory** with any ASE calculator (including foundation machine learning potentials such as **GRACE**, **MACE**, **CHGNet**, analytical potentials like **EMT**, or DFT calculators like **GPAW**, **VASP**, and **Quantum ESPRESSO**).

---

## Key Features

- **Works with ANY ASE Calculator:** Simply attach `atoms.calc = calc`.
- **Zero VASP Dependency:** Runs entirely in Python/RAM without dumping dozens of temporary directories or POSCAR/INCAR files to disk.
- **CLI Included:** Run directly from your command line: `elastic-ase POSCAR -c grace -m GRACE-2L-OAM`.
- **Automatic Symmetry Detection:** Automatically determines space group and crystal symmetry using `spglib`.
- **Supports All 7 Crystal Systems & Precise Sub-Classes:**
  - `cubic` (3 constants: $C_{11}, C_{12}, C_{44}$)
  - `hexagonal` (5 constants: $C_{11}, C_{12}, C_{13}, C_{33}, C_{44}$)
  - `tetragonal_1` (Laue $4/mmm$, 6 constants: $C_{11}, C_{12}, C_{13}, C_{33}, C_{44}, C_{66}$)
  - `tetragonal_2` (Laue $4/m$, 7 constants: including $C_{16}$)
  - `orthogonal` / `orthorhombic` (9 constants)
  - `trigonal_1` & `trigonal_2` (6 / 7 constants, handling $31m$ settings and arbitrary in-plane angles)
  - `monoclinic` (13 constants)
  - `triclinic` (21 constants)
- **Fast Symmetry & Orientation Check:** `elastic-ase --check POSCAR` provides instant space group classification, independent constants, and recommended workflow without running any energy evaluations.
- **Publication-Quality Anisotropic Visualizations:** Seamlessly integrated with `elastic-vis` for 2D polar plots and 3D anisotropic surfaces (Young's modulus, Shear modulus, Poisson's ratio, and linear compressibility).
- **Direct Stiffness Export:** Export plain $6 \times 6$ Voigt stiffness matrix files (`-o matrix.txt`), fully compatible with `elastic-vis` and `numpy.loadtxt`.
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
- `-c, --calculator`: `emt` or `grace` (default: `emt`).
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

---

## Crystallographic Conventions & Common Pitfalls

To prevent silent errors common in legacy elastic calculation codes, `ElasticASE` implements automatic validation and robust symmetry routing:

1. **Trigonal Systems & "31m"-Type Groups ($P\bar{3}1m, P\bar{3}1c, P31m, P31c, P312, P3_1 12, P3_2 12$):**
   - In standard settings with $a \parallel x$, the 2-fold axes in 31m-type groups are perpendicular to $a$, causing the stiffness matrix to have $C_{15}$ rather than $C_{14}$.
   - **Solution in ElasticASE:** Automatic symmetry detection routes all trigonal crystals to `trigonal_2`. This simultaneously fits both $C_{14}$ and $C_{15}$, making it completely invariant to in-plane rotation about $z$ with only 1 extra deformation mode.

2. **Rhombohedral Cells with 3-Fold Axis along $[111]$:**
   - Both `trigonal_1` and `trigonal_2` strain tensors require the 3-fold axis to lie along Cartesian $z$.
   - Primitive rhombohedral cells with equal angles (e.g., $a=b=c, \alpha=\beta=\gamma \ne 90^\circ$) have the 3-fold axis along $[111]$.
   - **Solution:** Convert the structure to the conventional hexagonal setting ($c \parallel z$) using `spglib` or ASE before running `ElasticASE`. A warning is automatically raised if an unconverted rhombohedral cell is detected.

3. **Monoclinic Structures (Unique Axis $b$ vs $c$):**
   - `ElasticASE` assumes standard International Tables convention: **unique axis $b$** ($\alpha = \gamma = 90^\circ$, $\beta \ne 90^\circ$).
   - If a structure is defined with unique axis $c$ ($\gamma \ne 90^\circ$), `ElasticASE` issues a warning to permute the axes before running.

4. **Hexagonal or Trigonal Cells with $c$ along $x$:**
   - Hexagonal and trigonal deformation matrices require the principal rotation axis $c$ to be parallel to Cartesian $z$ (`cell[2] || z`).
   - `ElasticASE` automatically checks the projection of the $c$-vector on $z$ and issues a warning if misalignment is detected.

