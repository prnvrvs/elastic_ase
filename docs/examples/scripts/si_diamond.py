"""Elastic constants of silicon (diamond) (Laue class m-3m) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python si_diamond.py
"""
from ase.build import bulk
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = bulk("Si", "diamond", a=5.43)
check_symmetry(atoms)

# Step 2: relax
atoms.calc = grace_fm("GRACE-2L-OAM")
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-4)

# Step 3: elastic constants
ec = ElasticConstants(atoms, fmax=1e-4)
results = ec.calculate()

# Step 4: check and save
ec.plot_fits("si_diamond_fits.png")
ec.save("si_diamond.json")
