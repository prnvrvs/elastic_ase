"""Elastic constants of a model crystal in P6/m (Laue group 6/m) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python model_p6m.py
"""
from ase.constraints import FixAtoms, FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["Ca", "Si", "O"], [(0, 0, 0), (1/3, 2/3, 0.5), (0.37, 0.12, 0.23)],
                spacegroup=175, cellpar=[6.4, 6.4, 3.9, 90, 90, 120])
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
ec.plot_fits("model_p6m_fits.png")
ec.save("model_p6m.json")
