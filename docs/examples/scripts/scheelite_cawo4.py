"""Elastic constants of scheelite CaWO4 (Laue group 4/m) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python scheelite_cawo4.py
"""
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["Ca", "W", "O"], [(0, 0.25, 0.625), (0, 0.25, 0.125), (0.15, 0.0075, 0.2106)],
                spacegroup=88, setting=2, cellpar=[5.243, 5.243, 11.376, 90, 90, 90],
                primitive_cell=True)
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
ec.plot_fits("scheelite_cawo4_fits.png")
ec.save("scheelite_cawo4.json")
