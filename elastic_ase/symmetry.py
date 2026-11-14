"""Symmetry classification and reporting for crystal structures.

Answers questions about a crystal before running calculations:
- Space group number and symbol
- Crystal system and Laue class
- Number of independent elastic constants in this orientation
- Which symmetry workflow is recommended (and whether the orientation matches)
- Works with both file paths (CIF, POSCAR, etc.) and `ase.Atoms` objects.
"""
from dataclasses import dataclass, field

import ase.io
import numpy as np
from ase import Atoms

from elastic_ase import auto, orientation, _spglib
from elastic_ase.systems import SYSTEMS

# (last space group, crystal system, Laue class, named workflow)
_RANGES = [
    (2, "triclinic", "-1", "triclinic"),
    (15, "monoclinic", "2/m", "monoclinic"),
    (74, "orthorhombic", "mmm", "orthorhombic"),
    (88, "tetragonal II", "4/m", "tetragonal_2"),
    (142, "tetragonal I", "4/mmm", "tetragonal_1"),
    (148, "trigonal II", "-3", "trigonal_2"),
    (167, "trigonal I", "-3m", "trigonal_1"),
    (194, "hexagonal", "6/mmm", "hexagonal"),
    (230, "cubic", "m-3m", "cubic"),
]

# trigonal I groups whose 2-fold axes are perpendicular to a in standard setting
# ("31m" type): C15 instead of C14 appears, so trigonal_2 is the named choice
TRIGONAL_31M = frozenset({149, 151, 153, 157, 159, 162, 163})


def classify(number: int):
    """(crystal_system, laue_class, named_workflow) for space group 1-230."""
    if not 1 <= number <= 230:
        raise ValueError(f"Space group number must be 1-230, got {number}")
    for last, system, laue, named in _RANGES:
        if number <= last:
            if named == "trigonal_1" and number in TRIGONAL_31M:
                named = "trigonal_2"
            return system, laue, named


@dataclass
class SymmetryReport:
    """Result of check_symmetry()."""

    source: str
    formula: str
    natoms: int
    space_group: str
    number: int
    crystal_system: str
    laue_class: str
    constants: tuple
    recommended: str
    recommended_works: bool
    working: list = field(default_factory=list)
    symprec: float = auto.SYMPREC

    @property
    def n_constants(self):
        return len(self.constants)

    def __getitem__(self, key):
        return getattr(self, key)

    def summary(self):
        lines = [
            f"{self.source}: {self.formula} ({self.natoms} atoms)",
            f"Space group:        {self.space_group} – {self.crystal_system} (Laue class {self.laue_class}), spglib symprec = {self.symprec} A",
            f"Elastic constants:  {self.n_constants} independent ({', '.join(self.constants)} in this orientation)",
        ]
        if self.recommended_works:
            lines.append(f"Recommended:        symmetry='{self.recommended}'  – works in this orientation")
            others = [w for w in self.working if w != self.recommended]
            if others:
                lines.append(f"Also works:         {', '.join(others)} (more calculations)")
        else:
            lines.append(
                f"Recommended:        symmetry='{self.recommended}'  – but NOT in this orientation: "
                "crystal axes are not aligned with standard Cartesian conventions."
            )
            if self.working:
                lines.append(f"Works as it is:     {', '.join(self.working)}")
            lines.append(
                f"To use symmetry='{self.recommended}', align the structure first or use symmetry='auto'."
            )
        lines.append(f"Any orientation:    symmetry='auto' – {self.n_constants} deformations")
        return "\n".join(lines)


def check_symmetry(structure, symprec=auto.SYMPREC, verbose=True):
    """Report crystal symmetry and recommended workflow for an Atoms object or structure file.

    Parameters
    ----------
    structure : ase.Atoms or str or os.PathLike
        ASE Atoms object or path to structure file (POSCAR, CIF, etc.).
    symprec : float
        spglib tolerance in Angstroms.
    verbose : bool
        If True, prints formatted summary.

    Returns
    -------
    SymmetryReport
    """
    if isinstance(structure, Atoms):
        atoms = structure
        source_name = f"Atoms({atoms.get_chemical_formula()})"
    else:
        source_name = str(structure)
        atoms = ase.io.read(source_name)

    spg_cell = (np.asarray(atoms.get_cell(), float), atoms.get_scaled_positions(), atoms.numbers)
    dataset = _spglib.call("get_symmetry_dataset", spg_cell, symprec=symprec)
    if dataset is None:
        raise ValueError(f"spglib could not determine symmetry for {source_name} (symprec={symprec})")

    number = dataset.number if hasattr(dataset, "number") else dataset["number"]
    system, laue, named = classify(number)
    rotations = auto.crystal_rotations(atoms, symprec)
    basis = auto.invariant_basis(rotations)

    working = sorted(
        (name for name, table in SYSTEMS.items() if orientation.table_covers(table, basis)),
        key=lambda name: (len(SYSTEMS[name].constants), name)
    )

    report = SymmetryReport(
        source=source_name,
        formula=atoms.get_chemical_formula(),
        natoms=len(atoms),
        space_group=orientation.space_group(atoms, symprec),
        number=int(number),
        crystal_system=system,
        laue_class=laue,
        constants=tuple(basis),
        recommended=named,
        recommended_works=named in working,
        working=working,
        symprec=symprec,
    )

    if verbose:
        print(report.summary())
    return report
