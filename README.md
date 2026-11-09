# ElasticASE

**ElasticASE** is a lightweight Python package for calculating second-order elastic constants ($C_{ij}$) using the **Atomic Simulation Environment (ASE)**.

It is the direct ASE counterpart and modern successor to `elastic_vasp`, operating **100% in-memory** with any ASE calculator (including machine learning potentials such as **GRACE**, **MACE**, **PACE**, **GAP**, **CHGNet**, or DFT calculators like **GPAW**, **VASP**, **Quantum ESPRESSO**).

---

## Key Features
- **Works with ANY ASE Calculator:** `atoms.calc = calc`
- **Zero VASP Dependency:** Runs entirely in Python/RAM without generating dozens of POSCAR/INCAR files or folders.
- **Supports All 7 Crystal Symmetries:**
  - `cubic` (3 constants: $C_{11}, C_{12}, C_{44}$)
  - `hexagonal` (5 constants: $C_{11}, C_{12}, C_{13}, C_{33}, C_{44}$)
  - `orthogonal` (orthorhombic / tetragonal, up to 9 constants)
  - `trigonal_1` & `trigonal_2` (6 / 7 constants)
  - `monoclinic` (13 constants)
  - `triclinic` (21 constants)
- **Automatic Symmetry Detection:** Via `spglib`.
- **Born Mechanical Stability Validation:** Evaluates stability criteria automatically.
- **Polycrystalline Averages:** Computes Voigt, Reuss, and Hill Bulk ($K$), Shear ($G$), Young's ($E$) moduli, Poisson's ratio ($\nu$), Pugh ratio ($K/G$), Cauchy pressure, and Zener anisotropy.
- **Optional Internal Atomic Relaxation:** Supports an ASE optimizer (e.g. `optimizer=BFGS`) for non-affine multi-atom crystal bases.

---

## Installation
```bash
cd elastic_ase
pip install -e .
```

---

## Quick Example (with GRACE)

```python
from ase.build import bulk
from tensorpotential.calculator import grace_fm
from elastic_ase import ElasticASE

# 1. Setup structure and ASE calculator
calc = grace_fm("GRACE-2L-OAM")
al = bulk("Al", "fcc", a=4.0428, cubic=True)
al.calc = calc

# 2. Run ElasticASE
el = ElasticASE(al, symmetry="cubic")
results = el.calculate()

print("C11:", results["stiffness_matrix"][0, 0])
print("C12:", results["stiffness_matrix"][0, 1])
print("C44:", results["stiffness_matrix"][3, 3])
```
