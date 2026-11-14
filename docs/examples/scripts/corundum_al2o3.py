"""Elastic constants of corundum Al2O3 (Laue group -3m) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python corundum_al2o3.py
"""
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["Al", "O"], [(0, 0, 0.352), (0.306, 0, 0.25)], spacegroup=167,
                cellpar=[4.759, 4.759, 12.99, 90, 90, 120], primitive_cell=True)
check_symmetry(atoms)

# Step 2: relax
atoms.calc = grace_fm("GRACE-2L-OAM")
atoms.set_constraint(FixSymmetry(atoms))
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-4)
atoms.set_constraint()

# Step 3: elastic constants
ec = ElasticConstants(atoms, fmax=1e-4)
results = ec.calculate()

# Step 4: check and save
ec.plot_fits("corundum_al2o3_fits.png")
ec.save("corundum_al2o3.json")
