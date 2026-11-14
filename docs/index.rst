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
       <a class="ev-btn" href="examples/step_by_step.html">Examples</a>
     </div>
   </section>

Two ways to run
---------------

.. list-table::
   :widths: 30 70
   :header-rows: 0

   * - **Any ASE calculator**
       (EMT, EAM, machine-learning potentials, GPAW, …)
     - ``ElasticConstants(atoms).calculate()`` – everything in memory, in one call.
       See :doc:`quickstart` and :doc:`guide/calculators`.
   * - **DFT run folders**
       (VASP, Quantum ESPRESSO)
     - ``ElasticConstants(atoms, code="vasp")`` writes one folder per calculation, runs them or
       lets you submit them, tracks which are done / failed / pending and reads the outputs
       back. See :doc:`guide/dft`.

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

Where to go
-----------

* **New here?** :doc:`installation`, the :doc:`quickstart`, then a complete calculation
  :doc:`step by step <examples/step_by_step>`.
* **Your crystal?** :doc:`examples/laue` – one worked example for each of the 11 Laue groups:
  the method elastic_ase uses and the result.
* **Your calculator?** :doc:`examples/calculators` (EAM, GRACE, MACE, GPAW, Quantum ESPRESSO)
  and :doc:`guide/dft` for VASP and Quantum ESPRESSO run folders.
* **Before production runs:** :doc:`guide/convergence`.
* **Something specific?** The :doc:`command line <guide/cli>`, the
  :doc:`Python API <api/index>`, or the :doc:`theory`.
* **How accurate is it?** :doc:`testing`. **Using it in a paper?** :doc:`citing`.

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

   guide/convergence
   guide/calculators
   guide/dft
   guide/orientation
   guide/results
   guide/cli

.. toctree::
   :maxdepth: 1
   :caption: Examples
   :hidden:

   examples/step_by_step
   examples/laue
   examples/calculators

.. toctree::
   :maxdepth: 1
   :caption: Reference
   :hidden:

   theory
   guide/systems
   api/index
   testing
   citing
   changelog
