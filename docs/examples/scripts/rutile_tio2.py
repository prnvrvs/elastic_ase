"""Elastic constants of rutile TiO2 (Laue class 4/mmm) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python rutile_tio2.py
"""
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["Ti", "O"], [(0, 0, 0), (0.305, 0.305, 0)], spacegroup=136,
                cellpar=[4.594, 4.594, 2.959, 90, 90, 90])
check_symmetry(atoms)

# Step 2: relax
atoms.calc = grace_fm("GRACE-2L-OAM")
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-4)

# Step 3: elastic constants
ec = ElasticConstants(atoms, fmax=1e-4)
results = ec.calculate()

# Step 4: check and save
ec.plot_fits("rutile_tio2_fits.png")
ec.save("rutile_tio2.json")
