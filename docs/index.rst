elastic_ase
===========

.. raw:: html

   <section class="ev-hero">
     <div class="ev-hero__eyebrow">elastic_ase documentation</div>
     <h1>Elastic constants with ASE calculators.</h1>
     <p>
       Give it a relaxed structure. elastic_ase finds the symmetry, applies the strains,
       gets the energies &ndash; from an ASE calculator, or from VASP and Quantum ESPRESSO
       run folders &ndash; and returns the 6&times;6 stiffness matrix.
     </p>
     <div class="ev-hero__actions">
       <a class="ev-btn ev-btn--primary" href="quickstart.html">Start here</a>
       <a class="ev-btn" href="installation.html">Install</a>
       <a class="ev-btn" href="guide/dft.html">VASP / Quantum ESPRESSO</a>
     </div>
   </section>

Two ways to run
---------------

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Fast calculators** (tested with EMT, EAM, GRACE and MACE; other ML potentials
       should work the same way)
     - ``ElasticConstants(atoms).calculate()`` – everything in memory.
       See :doc:`quickstart` and :doc:`guide/calculators`.
   * - **DFT codes** (VASP; Quantum ESPRESSO not yet tested with a real run)
     - ``ElasticConstants(atoms, code="vasp")`` writes one folder per calculation, runs or lets
       you submit them, tracks which are done / failed / pending and reads the outputs back.
       See :doc:`guide/dft`.

Both give the same results: stiffness matrix, independent constants, Born stability,
Voigt–Reuss–Hill averages, sound velocities, Debye temperature and more
(:doc:`guide/results`).

What elastic_ase checks for you
-------------------------------

* **Symmetry and orientation** – the crystal system is detected with spglib and the structure
  is rotated rigidly into the standard setting (your cell, atoms, magnetic moments and k-point
  mesh stay as they are). Without rotation an exact check warns when a deformation table
  cannot measure the crystal in its orientation (:doc:`guide/orientation`).
* **Residual stress** – a reference cell that is not relaxed biases the constants; the residual
  pressure is estimated from the fits and reported.
* **Every run** – crashed, unconverged (VASP, Quantum ESPRESSO), unreadable or mixed-up
  outputs are listed and can be repeated alone.
* **The curves** – ``plot_fits()`` shows every energy–strain curve with its fit.

What has and has not been tested is listed in :doc:`testing`.

Where to go
-----------

* **New here?** :doc:`installation`, then the :doc:`quickstart`.
* **DFT?** :doc:`guide/dft` – VASP and Quantum ESPRESSO, on your machine or on a cluster.
* **Something specific?** The :doc:`command line <guide/cli>`, the
  :doc:`Python API <api/index>`, or the :doc:`theory`.
* **Using it in a paper?** See :doc:`citing`.

elastic_ase is the ASE counterpart of `elastic_vasp <https://prnvrvs.github.io/elastic_vasp/>`_:
the same deformation tables and analysis, with energies from an ASE calculator or DFT run folders.

.. toctree::
   :maxdepth: 1
   :caption: Getting started
   :hidden:

   installation
   quickstart

.. toctree::
   :maxdepth: 1
   :caption: User guide
   :hidden:

   guide/calculators
   guide/dft
   guide/convergence
   guide/orientation
   guide/systems
   guide/results
   guide/cli

.. toctree::
   :maxdepth: 1
   :caption: Reference
   :hidden:

   theory
   testing
   citing
   api/index
   changelog
