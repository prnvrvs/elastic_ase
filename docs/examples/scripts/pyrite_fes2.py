"""Elastic constants of pyrite FeS2 (Laue group m-3) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python pyrite_fes2.py
"""
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["Fe", "S"], [(0, 0, 0), (0.385, 0.385, 0.385)], spacegroup=205,
                cellpar=[5.418, 5.418, 5.418, 90, 90, 90])
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
ec.plot_fits("pyrite_fes2_fits.png")
ec.save("pyrite_fes2.json")
