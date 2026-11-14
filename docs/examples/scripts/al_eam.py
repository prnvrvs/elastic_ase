"""Elastic constants of fcc aluminium with EAM potential of Mishin et al. (1999).

Needs: pip install elastic-ase, and the file Al99.eam.alloy from
       https://www.ctcms.nist.gov/potentials/ (aluminium, Mishin et al. 1999)
Run:   python al_eam.py
"""
from ase.build import bulk
from ase.calculators.eam import EAM
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS

from elastic_ase import ElasticConstants

calc = EAM(potential="Al99.eam.alloy")

# Structure: 1-atom fcc cell
atoms = bulk("Al", "fcc", a=4.04)

# Step 1: relax
atoms.calc = calc
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-4)

# Step 2: elastic constants
ec = ElasticConstants(atoms)
results = ec.calculate()
ec.save("al_eam.json")
