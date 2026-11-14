Python API
==========

The main entry point is ``elastic_ase.ElasticConstants``, the same class as
:class:`elastic_ase.ElasticASE <elastic_ase.calculator.ElasticASE>`; every analysis ends in an
:class:`~elastic_ase.results.ElasticResults`. The class delegates to small modules that can also
be used on their own:

* :mod:`elastic_ase.analysis` – from energies to results (pure functions),
* :mod:`elastic_ase.rundir` – run folders, saved energies, reading outputs,
* :mod:`elastic_ase.codes` – one adapter per DFT code (VASP, QE),
* :mod:`elastic_ase.systems` – the deformation tables,
* :mod:`elastic_ase.auto`, :mod:`elastic_ase.orientation`, :mod:`elastic_ase.symmetry` –
  symmetry, orientation and the automatic mode,
* :mod:`elastic_ase.average`, :mod:`elastic_ase.properties` – averages and derived properties,
* :mod:`elastic_ase.report`, :mod:`elastic_ase.plotting` – printed report and fit plots.

.. toctree::
   :maxdepth: 1

   calculator
   results
   analysis
   rundir
   codes
   properties
   average
   symmetry
   systems
   auto
   orientation
   solvers
   report
   plotting
