"""Elastic constants of baddeleyite ZrO2 (Laue group 2/m) with the GRACE-2L-OAM potential.

Needs: pip install elastic-ase tensorpotential
Run:   python baddeleyite_zro2.py
"""
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS
from ase.spacegroup import crystal
from tensorpotential.calculator import grace_fm

from elastic_ase import ElasticConstants, check_symmetry

# Step 1: structure
atoms = crystal(["Zr", "O", "O"], [(0.2758, 0.0411, 0.2082), (0.0703, 0.3359, 0.3406),
                                   (0.4423, 0.7549, 0.4789)], spacegroup=14,
                cellpar=[5.1505, 5.2116, 5.3173, 90, 99.23, 90])
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
ec.plot_fits("baddeleyite_zro2_fits.png")
ec.save("baddeleyite_zro2.json")
