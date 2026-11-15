Fast calculators (in memory)
============================

With any ASE calculator – EMT, a classical potential, a machine-learning potential, a DFT
code driven through ASE – everything runs in one call:

.. code-block:: python

   from elastic_ase import ElasticConstants

   atoms.calc = calc                  # any ASE calculator
   ec = ElasticConstants(atoms)
   results = ec.calculate()

The defaults – automatic symmetry, strains from −1 % to +1 % – can be changed with the
options below.

Before you start
----------------

**Converge the settings** of your calculator on the elastic constants themselves (k-mesh,
smearing, cutoff; see :doc:`convergence`), then **relax the structure** – cell and atoms –
with the same calculator and settings. The method assumes a stress-free reference: EMT
aluminium with 1.6 GPa of residual tension (lattice constant 1.4 % too large) gives
:math:`C_{11} = 47.6` instead of 53.2 GPa and :math:`C_{44} = 31.0` instead of 36.2 GPa.
elastic_ase estimates the residual stress from the fits and warns when it exceeds 0.5 GPa
along any deformation (see :doc:`results`).

.. code-block:: python

   from ase.filters import FrechetCellFilter
   from ase.optimize import BFGS

   BFGS(FrechetCellFilter(atoms)).run(fmax=1e-4)

Worked examples: :doc:`../examples/step_by_step`, one for every crystal symmetry in
:doc:`../examples/laue`, and one for each calculator in :doc:`../examples/calculators`.

Options
-------

All options are arguments of ``ElasticConstants``; leave out the ones you do not change.
Below is every argument with its default value, grouped by what it controls:

.. code-block:: python

   from ase.optimize import BFGS
   from elastic_ase import ElasticConstants

   ec = ElasticConstants(
       atoms,                       # relaxed structure, with atoms.calc attached

       # which deformations are applied
       symmetry="auto",             # symmetry table: detected, or e.g. "cubic"
       strains=[-0.01, -0.005, 0, 0.005, 0.01],  # strain values of every deformation

       # what happens inside every strained cell
       relaxation=True,             # relax the atoms (relaxed-ion constants)
       optimizer=BFGS,              # optimizer of that relaxation
       fmax=1e-3,                   # its force criterion (eV/Å)

       # orientation and symmetry detection
       reorient=True,               # rotate the structure into the standard setting
       frame="cartesian",           # frame of results["stiffness_matrix"]
       symprec=1e-3,                # spglib tolerance (Å)

       # output on the screen
       verbose=True,                # print progress and the report

       # DFT run folders only (see the DFT codes page)
       code=None,                   # "VASP" or "QE"
       command=None,                # e.g. "mpirun -np 16 vasp_std"
       output=None,                 # output file to read, e.g. "vasprun.xml"
   )
   results = ec.calculate(workdir=None)   # workdir="runs": keep every run, restartable

Two typical changes:

.. code-block:: python

   # accurate relaxation of the atoms in a multi-atom cell
   ec = ElasticConstants(atoms, fmax=1e-4)

   # clamped-ion constants with larger strains
   ec = ElasticConstants(atoms, relaxation=False, strains=[-0.02, -0.01, 0, 0.01, 0.02])

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Argument
     - Meaning
   * - ``atoms``
     - The relaxed structure (``ase.Atoms``) with the calculator attached as
       ``atoms.calc``. It is copied; your ``atoms`` object is not changed.
   * - ``symmetry``
     - ``"auto"`` (default): detected with spglib, the fewest deformations for the crystal.
       Or a named table: ``"cubic"``, ``"hexagonal"``, ``"tetragonal_1"``,
       ``"tetragonal_2"``, ``"trigonal_1"``, ``"trigonal_2"``, ``"orthogonal"``
       (= ``"orthorhombic"``), ``"monoclinic"``, ``"triclinic"`` (see :doc:`systems`).
       ``check_symmetry(atoms)`` shows which tables work for your structure.
   * - ``strains``
     - Strain values :math:`\delta` applied with every deformation (at least 3 different
       ones; default ``[-0.01, -0.005, 0, 0.005, 0.01]``). Use values symmetric around zero;
       ±0.5–1 % is a good range for most solids. Smaller strains make the fit sensitive to
       numerical noise, larger ones bring in anharmonic terms (see
       :ref:`choosing-the-strain-range`).
   * - ``relaxation``
     - ``True`` (default): relax the atoms inside every strained cell – *relaxed-ion*
       constants, the ones to compare with experiment. The atoms of a multi-atom cell move
       when the cell is strained, also in diamond Si or hcp metals under shear.
       ``False``: keep them at their scaled positions – *clamped-ion* constants.

       * The relaxation uses forces, so the calculator must provide them (energy-only
         calculators need ``relaxation=False``).
       * To keep some atoms in place, fix them with ``atoms.set_constraint(FixAtoms(...))``
         before creating ``ElasticConstants``; the relaxation leaves them where they are.
   * - ``optimizer``
     - The ASE optimizer class used for that relaxation (default ``BFGS``), e.g. ``FIRE`` or
       ``LBFGS`` from ``ase.optimize``. Pass the class, not an instance.
   * - ``fmax``
     - Force criterion of the relaxation in eV/Å (default ``1e-3``). ``1e-4`` gives more
       accurate constants for multi-atom cells at the cost of more force calls.
   * - ``reorient``
     - ``True`` (default): rotate the structure rigidly into the standard crystallographic
       setting before straining (cell shape and atoms unchanged). ``False``: keep it as
       given (see :doc:`orientation`).
   * - ``frame``
     - Frame of ``results["stiffness_matrix"]``: ``"cartesian"`` (default) – the frame of
       your input structure; ``"standard"`` – the standard crystallographic setting, the one
       tables of elastic constants use (see :doc:`orientation`).
   * - ``symprec``
     - Tolerance of the symmetry detection in Å (default ``1e-3``). Increase it (e.g.
       ``1e-2``) for a structure that is only approximately symmetric, such as one read from
       a file with few digits.
   * - ``verbose``
     - ``True`` (default): print progress and the full report; ``False``: no output on the
       screen. The messages also go to the ``"elastic_ase"`` logger.
   * - ``code``
     - DFT run folders only: ``"VASP"`` or ``"QE"``. Sets what
       ``write_inputs`` writes, which output file is read and checked, and the default
       command of ``run_jobs`` (see :doc:`dft`).
   * - ``command``
     - DFT run folders only: the shell command ``run_jobs`` runs in every folder, e.g.
       ``"mpirun -np 16 vasp_std"``.
   * - ``output``
     - DFT run folders only: the output file to read in every folder when it is not the
       code's default, e.g. ``"vasprun.xml"`` (any file ``ase.io.read`` reads).
   * - ``calculate(workdir=...)``
     - ``None`` (default): everything in memory. A folder name such as ``"runs"``: every
       energy is saved as it is computed and a repeated call only does what is missing (see
       `Keeping the runs: workdir`_ below).

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
  ``runs/<deformation>/strain_<value>/`` (for example ASE's ``Vasp`` or ``Espresso``);
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

Saved energies are only reused for the same run: each entry stores the strained cell, the
atoms (elements and positions before relaxation) and the relaxation settings
(``relaxation``, ``optimizer``, ``fmax``), and is ignored if any of them differ. The
calculator and its settings are not stored: use a new ``workdir`` for every calculator
setting, for example in a convergence test.

A typical loop for a calculator that sometimes fails:

.. code-block:: python

   for attempt in range(3):
       if ec.is_complete("runs"):
           break
       ec.run_pending("runs")
   results = ec.calculate(workdir="runs")
