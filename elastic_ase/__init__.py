__version__ = "0.0.1"

from elastic_ase.calculator import ElasticASE
from elastic_ase.properties import derived_properties
from elastic_ase.results import ElasticResults

from elastic_ase.symmetry import check_symmetry, classify, SymmetryReport
from elastic_ase.deformations import get_deformation_specs
from elastic_ase.solvers import solve_stiffness, check_mechanical_stability
from elastic_ase.average import calculate_polycrystalline_averages
from elastic_ase.systems import SYSTEMS, CrystalSystem
from elastic_ase import auto

# Same class under a descriptive name: ec = ElasticConstants(atoms)
ElasticConstants = ElasticASE

__all__ = [
    "ElasticASE",
    "ElasticConstants",
    "ElasticResults",
    "derived_properties",
    "check_symmetry",
    "classify",
    "SymmetryReport",
    "get_deformation_specs",
    "solve_stiffness",
    "check_mechanical_stability",
    "calculate_polycrystalline_averages",
    "SYSTEMS",
    "CrystalSystem",
    "auto",
]
