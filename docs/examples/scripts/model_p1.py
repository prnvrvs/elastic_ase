"""Elastic constants of a model crystal in P-1 (Laue group -1) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python model_p1.py
"""
from ase.constraints import FixAtoms, FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["Au", "Cu"], [(0, 0, 0), (0.13, 0.27, 0.41)], spacegroup=2,
                cellpar=[3.9, 4.3, 4.7, 82, 97, 104])
check_symmetry(atoms)

# Step 2: relax
atoms.calc = grace_fm("GRACE-2L-OAM")
atoms.set_constraint([FixAtoms(indices=range(len(atoms))), FixSymmetry(atoms)])
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-3)
atoms.set_constraint()

# Step 3: elastic constants
ec = ElasticConstants(atoms, relaxation=False)  # model: atoms stay fixed
results = ec.calculate()

# Step 4: check and save
ec.plot_fits("model_p1_fits.png")
ec.save("model_p1.json")
