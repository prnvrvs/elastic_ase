"""Elastic constants of fcc aluminium with GPAW (DFT, PBE, plane waves).

Needs: pip install elastic-ase gpaw (with its PAW datasets)
Run:   mpiexec -n 6 gpaw python al_gpaw.py   (or: python al_gpaw.py)
"""
from ase.build import bulk
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from gpaw import GPAW, PW, MethfesselPaxton

from elastic_ase import ElasticConstants

calc = GPAW(mode=PW(400), xc="PBE", kpts={"size": (32, 32, 32), "gamma": True},
            occupations=MethfesselPaxton(0.2, order=1),
            convergence={"energy": 1e-5, "forces": 1e-3}, txt="al_gpaw.txt")

# Structure: 1-atom fcc cell
atoms = bulk("Al", "fcc", a=4.04)

# Step 1: relax
atoms.calc = calc
BFGS(FrechetCellFilter(atoms), logfile=None).run(fmax=1e-4)

# Step 2: elastic constants
ec = ElasticConstants(atoms)
results = ec.calculate()
ec.save("al_gpaw.json")
