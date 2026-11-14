"""Elastic constants of black phosphorus (Laue group mmm) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python black_p.py
"""
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["P"], [(0, 0.1017, 0.0806)], spacegroup=64,
                cellpar=[3.314, 10.48, 4.376, 90, 90, 90], primitive_cell=True)
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
ec.plot_fits("black_p_fits.png")
ec.save("black_p.json")
