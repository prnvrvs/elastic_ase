"""Elastic constants of ilmenite FeTiO3 (Laue class -3) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python ilmenite_fetio3.py
"""
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["Fe", "Ti", "O"], [(0, 0, 0.3554), (0, 0, 0.1464), (0.317, 0.020, 0.245)],
                spacegroup=148, cellpar=[5.088, 5.088, 14.09, 90, 90, 120],
                primitive_cell=True)
check_symmetry(atoms)

# Step 2: relax
atoms.calc = grace_fm("GRACE-2L-OAM")
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-4)

# Step 3: elastic constants
ec = ElasticConstants(atoms, fmax=1e-4)
results = ec.calculate()

# Step 4: check and save
ec.plot_fits("ilmenite_fetio3_fits.png")
ec.save("ilmenite_fetio3.json")
