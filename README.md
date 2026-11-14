# elastic_ase

Second-order elastic constants ($C_{ij}$) from energy–strain calculations with ASE
calculators. elastic_ase is the ASE counterpart of
[elastic_vasp](https://github.com/prnvrvs/elastic_vasp): the same deformation tables and
analysis, with the energies from an ASE calculator (run in memory) or from DFT run folders
(VASP, Quantum ESPRESSO) that elastic_ase writes and reads back.

**Documentation:** https://prnvrvs.github.io/elastic_ase/

## Features

- Symmetry detection with spglib and deformation tables for all crystal systems: cubic,
  hexagonal, tetragonal (4/mmm and 4/m), trigonal (with $C_{14}$ and $C_{15}$),
  orthorhombic, monoclinic, triclinic; an orientation-aware automatic mode.
- The structure is rotated rigidly into the standard setting (cell, atoms and per-atom data
  unchanged); without rotation, an exact check warns when a table does not fit the
  orientation.
- Energies from an ASE calculator, or run folders for a DFT code with done / failed /
  pending tracking and restarts.
- Results: stiffness matrix, Born stability, Voigt–Reuss–Hill averages, sound velocities,
  Debye temperature, residual pressure of the reference cell, energy–strain fit plots,
  save/load; directional properties and plots with
  [elastic-vis](https://github.com/prabhat7iitm/elastic-vis).
- Command line tool `elastic-ase`.

## Installation

```bash
pip install elastic-ase                                          # from PyPI
pip install "git+https://github.com/prnvrvs/elastic_ase.git"     # current version from GitHub
```

Python 3.10 or newer. Optional, for directional plots and extremes (not on PyPI):
`pip install "git+https://github.com/prabhat7iitm/elastic-vis.git"`.

## Quick start

With an ASE calculator (here EMT; relax the structure with the same calculator first):

```python
from ase.build import bulk
from ase.calculators.emt import EMT
from elastic_ase import ElasticConstants

atoms = bulk("Al", "fcc", a=3.9943, cubic=True)   # EMT equilibrium lattice constant
atoms.calc = EMT()
results = ElasticConstants(atoms).calculate()
print(results["independent_constants"])          # {'C11': ..., 'C12': ..., 'C44': ...}
```

With VASP run folders (run them with `run_jobs` or submit them on a cluster):

```python
from ase.io import read
from elastic_ase import ElasticConstants

ec = ElasticConstants(read("CONTCAR"), code="vasp", command="mpirun -np 16 vasp_std")
ec.write_inputs("runs", files={"INCAR": "INCAR", "KPOINTS": "KPOINTS", "POTCAR": "POTCAR"})
ec.run_jobs("runs")                   # only runs without a usable OUTCAR
ec.status("runs")                     # done / failed / pending
results = ec.read_results("runs")
```

From the command line: `elastic-ase POSCAR --check` (symmetry only) or
`elastic-ase POSCAR -c emt -o C.txt`.

The [documentation](https://prnvrvs.github.io/elastic_ase/) has a
[step-by-step example](https://prnvrvs.github.io/elastic_ase/examples/step_by_step.html),
[an example for every crystal symmetry](https://prnvrvs.github.io/elastic_ase/examples/laue.html)
and [for every calculator](https://prnvrvs.github.io/elastic_ase/examples/calculators.html),
guides to [VASP and Quantum ESPRESSO](https://prnvrvs.github.io/elastic_ase/guide/dft.html),
[convergence](https://prnvrvs.github.io/elastic_ase/guide/convergence.html),
[symmetry and orientation](https://prnvrvs.github.io/elastic_ase/guide/orientation.html) and
[the results](https://prnvrvs.github.io/elastic_ase/guide/results.html), the
[theory](https://prnvrvs.github.io/elastic_ase/theory.html), the
[accuracy](https://prnvrvs.github.io/elastic_ase/testing.html) against other methods and
experiment, and the [Python API](https://prnvrvs.github.io/elastic_ase/api/index.html).

Problems and questions: [issue tracker](https://github.com/prnvrvs/elastic_ase/issues).

## Development

```bash
git clone https://github.com/prnvrvs/elastic_ase.git
cd elastic_ase
pip install -e ".[test,doc]"
python -m unittest discover -s tests
make -C docs html
```

## Citing

If elastic_ase is useful for your work, please cite the papers listed in
[Citing](https://prnvrvs.github.io/elastic_ase/citing.html).

## License

MIT, see [LICENSE](https://github.com/prnvrvs/elastic_ase/blob/main/LICENSE).
