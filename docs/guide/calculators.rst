Fast calculators (in memory)
============================

With a calculator that returns an energy quickly – EMT, a classical potential, a
machine-learning potential – everything runs in one call (tested with EMT, GRACE and GPAW;
see :doc:`../testing`):

.. code-block:: python

   from elastic_ase import ElasticConstants

   atoms.calc = calc                  # any ASE calculator
   ec = ElasticConstants(atoms)
   results = ec.calculate()

The defaults – automatic symmetry, strains from −1 % to +1 % – can be changed with the
options below.

Before you start
----------------

**Relax the structure** – cell and atoms – with the same calculator and settings. The method
assumes a stress-free reference; residual stress changes the fitted curvatures by about its
own size. elastic_ase estimates the residual pressure from the fits and warns above 0.5 GPa
(see :doc:`results`).

.. code-block:: python

   from ase.filters import FrechetCellFilter
   from ase.optimize import BFGS

   BFGS(FrechetCellFilter(atoms)).run(fmax=1e-4)

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Argument
     - Meaning
   * - ``symmetry``
     - ``"auto"`` (default: detected with spglib) or a named table: ``"cubic"``,
       ``"hexagonal"``, ``"tetragonal_1"``, ``"tetragonal_2"``, ``"trigonal_1"``,
       ``"trigonal_2"``, ``"orthogonal"`` (= ``"orthorhombic"``), ``"monoclinic"``,
       ``"triclinic"`` (see :doc:`systems`).
   * - ``strains``
     - Strain values :math:`\delta` (at least 3 different ones; default
       ``[-0.01, -0.005, 0, 0.005, 0.01]``). Use values symmetric around zero; ±0.5–1 % is a
       good range for most solids, larger strains bring in anharmonic terms.
   * - ``optimizer``, ``fmax``
     - An ASE optimizer class (e.g. ``BFGS``) to relax the atoms inside every strained cell
       (needed when atoms are not fixed by symmetry, e.g. most multi-atom cells).
   * - ``reorient``
     - ``True`` (default): rotate the structure rigidly into the standard setting.
       ``False``: keep it as given (see :doc:`orientation`).
   * - ``frame``
     - ``"cartesian"`` (default): C in the frame of the calculation; ``"standard"``: in the
       standard crystallographic setting.
   * - ``symprec``
     - spglib tolerance in Å (default 1e-3).
   * - ``verbose``
     - Print the report (default ``True``).

Number of calculations
----------------------

One unstrained cell plus every deformation at every non-zero strain:
``len(ec.runs())``. For the default 5 strains that is 13 (cubic), 21 (hexagonal),
37 (orthorhombic) up to 85 (triclinic).

Keeping the runs: ``workdir``
-----------------------------

For slower calculators, or to keep every energy:

.. code-block:: python

   results = ec.calculate(workdir="runs")

* each energy is saved to ``runs/energies.json`` as soon as it is computed;
* a file-based calculator (one with a ``directory``) runs every calculation in its own folder
  ``runs/<deformation>/strain_<value>/`` (tested with ASE's ``Vasp`` running real VASP; not yet
  with ``Espresso``);
* a run that raises an error, or whose calculator reports ``converged = False``, is recorded
  in ``runs/failed.json``; the other runs still go ahead and a ``RuntimeError`` at the end
  lists the failures;
* calling ``calculate(workdir="runs")`` again repeats only what is missing.

The individual steps are also available. ``status`` returns
``{"done": [...], "failed": {run: reason}, "pending": [...]}``, ``missing`` the failed and
pending runs in run order, and ``run_pending`` runs them (or ``runs=[...]`` for chosen ones):

.. code-block:: python

   ec.status("runs")
   ec.missing("runs")
   ec.is_complete("runs")     # True when every run has an energy
   ec.run_pending("runs")
   ec.calculate(workdir="runs")

Saved energies are only reused for the same structure and strains: each entry stores its cell
and number of atoms and is ignored if they do not match.

A typical loop for a calculator that sometimes fails:

.. code-block:: python

   for attempt in range(3):
       if ec.is_complete("runs"):
           break
       ec.run_pending("runs")
   results = ec.calculate(workdir="runs")
