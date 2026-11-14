"""Elastic constants of fcc aluminium with the MACE-MP machine-learning potential (medium-mpa-0).

Needs: pip install elastic-ase mace-torch
Run:   python al_mace.py
"""
from ase.build import bulk
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from mace.calculators import mace_mp

from elastic_ase import ElasticConstants

calc = mace_mp(model="medium-mpa-0", default_dtype="float64")  # double precision

# Structure: 1-atom fcc cell
atoms = bulk("Al", "fcc", a=4.04)

# Step 1: relax
atoms.calc = calc
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-4)

# Step 2: elastic constants
ec = ElasticConstants(atoms)
results = ec.calculate()
ec.save("al_mace.json")
