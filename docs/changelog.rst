Changelog
=========

0.0.4
-----

* Fixed: with ``relaxation=True`` the unstrained cell is now relaxed too, to the same force
  criterion as the strained cells. Before, leftover forces in the input structure (e.g. from
  a relaxation with a looser ``fmax``) biased every constant – for a 4-atom Cu cell with
  forces up to 0.09 eV/Å, :math:`C_{44}` was 27 % too low.
* Fixed: energies saved in a ``workdir`` were reused for runs with other relaxation settings
  or other atom positions in the same cell. Each saved energy now also records the atoms and
  the relaxation settings (``relaxation``, ``optimizer``, ``fmax``) and is reused only if they
  match; energies saved by older versions count as unrelaxed runs.
* Fixed: under MPI every process wrote ``energies.json`` (and ``save``, ``export_matrix``,
  ``plot_fits``) at the same time, which could corrupt the file; now only rank 0 writes.
* Fixed: the command line rejected negative strains (``--strains -0.01,0,0.01``); a
  ``--plot`` run without elastic-vis and invalid strains now stop with a short error instead
  of a traceback. ``-o result.json`` now also holds the independent constants, derived
  properties, residual stress and every energy.
* Strain values that are equal to 6 decimals are rejected (they would share one run name).
* The residual-stress warning no longer says the bias is "a similar amount": the constants
  can be off by several times the residual stress.
* Requires ASE 3.23 or newer (the QE support uses ``EspressoProfile``); CI also tests Python 3.11.
* Messages name Quantum ESPRESSO "QE", e.g. "QE SCF not converged"; an unknown ``code=``
  lists ``['QE', 'VASP']``. The code name works in any case (``"VASP"``, ``"vasp"``, ``"QE"``,
  ``"qe"``); the documentation uses ``code="VASP"`` and ``code="QE"``.
* Documentation: every option of ``ElasticConstants`` in one annotated call
  (:doc:`guide/calculators`); how to choose the strain range, with a strain-range scan
  (:ref:`choosing-the-strain-range`); the matrix of every volume-conserving deformation
  (:doc:`theory`).

0.0.3
-----

* New: ``relaxation=True`` (the default) relaxes the atom positions inside every strained cell
  with BFGS, giving the relaxed-ion constants; ``relaxation=False`` keeps them fixed
  (clamped-ion). Before, the atoms were relaxed only when an ``optimizer`` was passed, so
  multi-atom cells gave clamped-ion constants by default. ``optimizer=`` still chooses the
  optimizer. Results of fcc and bcc metals do not change (symmetry keeps their atoms in place). The report shows the setting ("Atom relaxation");
  the command line has ``--no-relaxation`` and ``--fmax``.
* Fixed: an unconverged QE run whose output ASE cannot read was reported with
  ASE's read error; it is now reported as "SCF not converged" (the convergence check runs
  before reading).
* Fixed: a relative QE ``pseudo_dir`` (e.g. ``"pseudo"``) was looked up inside
  every run folder, so pw.x could not find the pseudopotentials; it is now taken relative to
  the folder you start from.
* Under MPI (e.g. GPAW with ``mpiexec``) the report is printed once, not by every process.
* ``independent_constants``, ``standard_constants`` and ``polycrystalline_averages`` hold plain
  Python floats, so printing them shows ``{'C11': 53.2, ...}`` instead of ``np.float64(...)``.
* Documentation reorganised around examples: a complete calculation step by step
  (:doc:`examples/step_by_step`), one example for each of the 11 Laue groups
  (:doc:`examples/laue`) and one for each calculator, including GPAW and QE
  (:doc:`examples/calculators`). The "What is tested" page is replaced by :doc:`testing`
  (comparison with the stress–strain method, published potentials and experiment).

0.0.2
-----

* Fixed: the Laue class reported for space groups 168–176 was "6/mmm" instead of "6/m", and for
  195–206 "m-3m" instead of "m-3" (``check_symmetry``, ``--check``). The deformation tables and
  results were not affected: both pairs share the same independent constants.
* Documentation: new :doc:`guide/convergence` page (what to converge, and a script for the
  k-mesh test); "What is tested" extended with all 11 Laue classes, 15 real materials compared
  with the stress–strain method, comparisons with experiment, MACE and GPAW.
* Command line: a ``-c module:function`` calculator that needs arguments gives a clear error
  instead of a traceback.
* Documentation: shorter README pointing to the docs, new :doc:`citing` page, links to the
  Python/NumPy/ASE documentation in the API reference, layout fixed for narrow screens, and
  the :doc:`testing` page extended (VASP, the Mishin 1999 Al EAM potential compared with its
  published constants, GRACE, GPAW).

0.0.1 (first release)
---------------------

* Run folders for DFT codes: ``code="VASP"`` / ``"QE"`` with ``write_inputs``, ``run_jobs``,
  ``status``, ``missing``, ``is_complete`` and ``read_results`` (plus ``vasp_*`` / ``qe_*``
  shortcuts); ``calculate(workdir=...)`` saves every energy and resumes after failures.
  Unconverged VASP and QE runs are detected.
* Rotation into the standard setting is now a rigid rotation of the given structure: the cell,
  the atoms and per-atom data (magnetic moments, tags, constraints) are kept. Before, the cell
  was replaced by spglib's conventional cell.
* Fixed: ``frame="standard"`` rotated the stiffness matrix twice for rotated input.
* Fixed: the stability check used named criteria in a rotated frame; it now uses the
  eigenvalues. Zener anisotropy and Cauchy pressure use the standard setting.
* Orientation check: exact test of whether a named table fits the crystal's orientation
  (replaces checks based on the cell shape).
* New: residual-pressure estimate from the fits, derived properties (density, sound
  velocities, Debye temperature, Kleinman parameter, hardness estimates), ``plot_fits``,
  ``save`` / ``load``, ``ElasticConstants`` alias.
* If spglib fails, the triclinic table is used instead of assuming cubic. At least three
  different strains are required.
* Restructured: ``calculator.py`` is split into :mod:`~elastic_ase.analysis`
  (energies to results), :mod:`~elastic_ase.rundir` (run folders), :mod:`~elastic_ase.codes`
  (one adapter per DFT code), :mod:`~elastic_ase.results` (``ElasticResults``),
  :mod:`~elastic_ase.report` and :mod:`~elastic_ase.plotting`; the public API is unchanged and
  all results are identical. The attributes ``ec.stiffness``, ``ec.averages``, … are now
  read-only views of ``ec.results``. ``deformations.py`` and ``solve_stiffness`` no longer keep
  their own copies of the tables. pandas is no longer needed. The command line accepts any
  calculator as ``-c module:function``.
* Packaging: a single ``pyproject.toml`` (MIT license) ready for PyPI; elastic-vis stays an
  optional separate install from GitHub (PyPI does not accept dependencies on GitHub URLs); Python 3.10 or newer; CI on Python 3.10, 3.12 and 3.13.

Before 0.0.1
------------

* Development versions: ASE port of elastic_vasp's deformation tables, symmetry checks, automatic
  mode, averages and command line tool.
