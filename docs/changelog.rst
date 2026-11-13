Changelog
=========

0.0.1 (first release)
---------------------

* Run folders for DFT codes: ``code="vasp"`` / ``"qe"`` with ``write_inputs``, ``run_jobs``,
  ``status``, ``missing``, ``is_complete`` and ``read_results`` (plus ``vasp_*`` / ``qe_*``
  shortcuts); ``calculate(workdir=...)`` saves every energy and resumes after failures.
  Unconverged VASP and Quantum ESPRESSO runs are detected.
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
* Restructured: ``calculator.py`` (1,429 lines) is split into :mod:`~elastic_ase.analysis`
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
