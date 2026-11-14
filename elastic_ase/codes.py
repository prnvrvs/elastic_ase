"""DFT code adapters: what to write into a run folder, which output to read, how to recognise an
unconverged run, and the default command.

Adding a code means adding one subclass of :class:`Code` and registering it in ``CODES``.
"""
import os
from dataclasses import dataclass, field


@dataclass
class WritePlan:
    """How the run folders are written."""
    use_calculator: bool                     # write the attached calculator's own input
    format: str = "vasp"                     # otherwise: ase.io.write format of the structure
    filename: str = "POSCAR"                 # ... and its file name
    notes: list = field(default_factory=list)  # messages for the user


class Code:
    """Base adapter. ``name``, ``input``, ``output`` and ``command`` describe the code."""
    name = ""
    input = ""
    output = ""
    command = ""

    def plan(self, atoms, calc, files):
        """WritePlan for these atoms, attached calculator and extra input files."""
        raise NotImplementedError

    def not_converged(self, path):
        """Reason string if the output at ``path`` belongs to this code and is not converged,
        else None."""
        return None


class Vasp(Code):
    name, input, output, command = "vasp", "POSCAR", "OUTCAR", "vasp_std"

    def plan(self, atoms, calc, files):
        from ase.calculators.vasp import Vasp as AseVasp
        if not files and isinstance(calc, AseVasp):
            return WritePlan(use_calculator=True)
        species = species_order(atoms)
        if "POTCAR" in files:
            elements = potcar_order(files["POTCAR"])
            if elements and elements != species:
                raise ValueError(f"POTCAR elements {elements} do not match the POSCAR species order {species}.")
        missing = [name for name in ("INCAR", "KPOINTS", "POTCAR") if name not in files]
        notes = [f"no {', '.join(missing)} given; add them to every folder before running VASP."] if missing else []
        return WritePlan(use_calculator=False, format="vasp", filename="POSCAR", notes=notes)

    def not_converged(self, path):
        if vasp_converged(path) is False:
            return "VASP electronic loop not converged (EDIFF not reached)"
        return None


class QuantumEspresso(Code):
    name, input, output = "qe", "espresso.pwi", "espresso.pwo"
    command = "pw.x -in espresso.pwi > espresso.pwo"

    def plan(self, atoms, calc, files):
        if calc is None or not hasattr(calc, "write_inputfiles"):
            raise ValueError(
                "QE input needs an ASE Espresso calculator attached (atoms.calc = "
                "Espresso(profile=..., pseudopotentials=..., input_data=..., kpts=...)): the "
                "pseudopotentials and settings come from it. For structures only, use "
                "write_inputs(workdir, format='espresso-in') without code."
            )
        return WritePlan(use_calculator=True)

    def not_converged(self, path):
        if qe_converged(path) is False:
            return "QE SCF not converged (convergence NOT achieved)"
        return None


CODES = {"vasp": Vasp(), "qe": QuantumEspresso()}
ALIASES = {"espresso": "qe", "quantum_espresso": "qe", "pwscf": "qe"}


def get_code(name):
    """The adapter for a code name ("VASP", "QE", "espresso", ...; any case), or None for None."""
    if name is None:
        return None
    if isinstance(name, Code):
        return name
    key = ALIASES.get(str(name).lower(), str(name).lower())
    if key not in CODES:
        raise ValueError(f"Unknown code {name!r}; choose from {sorted(c.upper() for c in CODES)} (or give output= and command=)")
    return CODES[key]


def not_converged(path, code=None):
    """Reason if the output is recognisably unconverged: by the given code's check, or, without
    a code, by any code that recognises the file. None otherwise."""
    adapters = [code] if code is not None else list(CODES.values())
    for adapter in adapters:
        reason = adapter.not_converged(path)
        if reason:
            return reason
    return None


# --- helpers ------------------------------------------------------------------------------

def species_order(atoms):
    """Element order of the POSCAR species line; error if a species is split into blocks."""
    order = []
    for symbol in atoms.get_chemical_symbols():
        if not order or order[-1] != symbol:
            order.append(symbol)
    if len(order) != len(set(order)):
        raise ValueError(
            f"The atoms are not grouped by element ({' '.join(order)}), so the POSCAR would list a "
            "species more than once and no POTCAR can match it. Sort them first, e.g. "
            "atoms = atoms[atoms.numbers.argsort(kind='stable')], then create ElasticConstants again."
        )
    return order


def potcar_order(path):
    """Element order of a (concatenated) POTCAR, from its TITEL lines (e.g. 'PAW_PBE Cu_pv ...')."""
    order = []
    with open(path, errors="replace") as f:
        for line in f:
            if line.strip().startswith("TITEL"):
                order.append(line.split("=", 1)[1].split()[1].split("_")[0])
    return order


def vasp_converged(path):
    """True/False for a VASP OUTCAR (electronic convergence), None for other files."""
    if not os.path.basename(path).startswith("OUTCAR"):
        return None
    try:
        from ase.calculators.vasp import Vasp as AseVasp
        with open(path, errors="replace") as f:
            return bool(AseVasp().read_convergence(lines=f.readlines()))
    except Exception:
        return None


def qe_converged(path):
    """True/False for a QE pw.x output (SCF convergence), None for other files."""
    try:
        with open(path, errors="replace") as f:
            text = f.read()
    except OSError:
        return None
    if "Program PWSCF" not in text:
        return None
    return "convergence NOT achieved" not in text
