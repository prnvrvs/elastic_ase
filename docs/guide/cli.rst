Command line
============

The ``elastic-ase`` command is a quick tool for fast calculators: EMT (for tests), GRACE
foundation models, or any ASE calculator given as ``module:function``. For DFT codes use
the Python API (:doc:`dft`), where the calculator settings, run folders and restarts are under
your control.

.. code-block:: bash

   elastic-ase POSCAR --check                # symmetry only
   elastic-ase POSCAR -c emt
   elastic-ase in.cif -c grace -m GRACE-2L-OAM -o results.json
   elastic-ase in.cif -c grace -o C.txt      # plain 6x6 matrix
   elastic-ase in.cif -c grace --plot --Young --Shear --format pdf

The structure can be any file ASE reads (POSCAR, CIF, extended XYZ, …).

Options
-------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Meaning
   * - ``--check``
     - Space group, independent constants and which tables fit; no calculation.
   * - ``-c``, ``--calculator``
     - ``emt`` (default), ``grace``, or ``module:function`` – a function (or class) called without
       arguments that returns an ASE calculator, e.g. ``ase.calculators.emt:EMT``. Tested with
       EMT only.
   * - ``-m``, ``--model``
     - GRACE model name (default ``GRACE-2L-OAM``).
   * - ``-s``, ``--symmetry``
     - ``auto`` (default) or a named table (:doc:`systems`).
   * - ``--strains``
     - Comma-separated strains (default ``-0.01,-0.005,0.0,0.005,0.01``).
   * - ``--frame``
     - ``cartesian`` (default) or ``standard`` (:doc:`orientation`).
   * - ``--no-reorient``
     - Keep the structure in its given orientation.
   * - ``-o``, ``--output``
     - Save the results: ``.json`` (everything), ``.csv`` (one row of moduli and
       :math:`C_{ij}`) or any other extension (plain 6×6 matrix).
   * - ``-q``, ``--quiet``
     - No report.

Plots (need elastic-vis): ``--plot`` with ``--Young``, ``--Shear``, ``--Poisson``, ``--LC``
(default Young and shear), ``--2D`` / ``--3D`` (default both), ``--format`` (png, pdf, svg),
``--output-prefix``, ``--style`` (JSON style file) and ``--show``.

``elastic-ase --help`` lists everything.
