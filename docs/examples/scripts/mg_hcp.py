"""Elastic constants of magnesium (hcp) (Laue group 6/mmm) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python mg_hcp.py
"""
from ase.build import bulk
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = bulk("Mg", "hcp", a=3.21, c=5.21)
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
ec.plot_fits("mg_hcp_fits.png")
ec.save("mg_hcp.json")
