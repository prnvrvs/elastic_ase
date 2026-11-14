Examples by calculator
======================

elastic_ase works with any ASE calculator. The examples on this page all compute fcc aluminium
(1-atom cell, cubic, 3 constants from 13 energies), so the calculators can be compared directly.
Each script has two steps:

* **Step 1: relax** – the structure, with the calculator you use for the elastic constants;
* **Step 2: elastic constants** – ``ElasticConstants(atoms).calculate()``.

Only the calculator line differs between the scripts. The full script for the EAM potential:

.. literalinclude:: scripts/al_eam.py
   :language: python
   :start-at: from ase.build import bulk


EAM (classical potential)
-------------------------

The aluminium potential of Mishin et al. (1999); download ``Al99.eam.alloy`` from the `NIST Interatomic Potentials Repository <https://www.ctcms.nist.gov/potentials/>`_. Its published constants are :math:`C_{11} = 114.0`, :math:`C_{12} = 61.6`, :math:`C_{44} = 31.6` GPa; with smaller strains (``strains=[-0.0025, -0.00125, 0, 0.00125, 0.0025]``) elastic_ase gives 113.5, 61.3 and 31.6 GPa – at the default ±1 % the anharmonic terms of this potential lower :math:`C_{11}` slightly.

**Screen output** (shortened):

.. literalinclude:: output/al_eam.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - 112.8
     - 60.5
     - 31.7

Script: :download:`al_eam.py <scripts/al_eam.py>`.

GRACE (machine-learning potential)
----------------------------------

A foundation model trained on PBE data; ``pip install tensorpotential``.

The calculator (the rest is as above):

.. literalinclude:: scripts/al_grace.py
   :language: python
   :start-after: from elastic_ase import ElasticConstants
   :end-before: # Structure

**Screen output** (shortened):

.. literalinclude:: output/al_grace.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - 112.8
     - 60.5
     - 33.9

Script: :download:`al_grace.py <scripts/al_grace.py>`.

MACE (machine-learning potential)
---------------------------------

``pip install mace-torch``. Use double precision: the energy differences between the strained cells are small.

The calculator (the rest is as above):

.. literalinclude:: scripts/al_mace.py
   :language: python
   :start-after: from elastic_ase import ElasticConstants
   :end-before: # Structure

**Screen output** (shortened):

.. literalinclude:: output/al_mace.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - 114.8
     - 64.9
     - 29.3

Script: :download:`al_mace.py <scripts/al_mace.py>`.

GPAW (DFT)
----------

PBE with plane waves. Metals need dense k-point meshes: for Al the shear constants still change by a few GPa between 30³ and 32³ k-points. Converge the mesh on the elastic constants (:doc:`../guide/convergence`). Run in parallel with ``mpiexec -n 6 python al_gpaw.py``; the report is printed once.

The calculator (the rest is as above):

.. literalinclude:: scripts/al_gpaw.py
   :language: python
   :start-after: from elastic_ase import ElasticConstants
   :end-before: # Structure

**Screen output** (shortened):

.. literalinclude:: output/al_gpaw.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - 118.8
     - 57.8
     - 33.9

Script: :download:`al_gpaw.py <scripts/al_gpaw.py>`.

Quantum ESPRESSO in run folders
-------------------------------

For DFT codes it is often better to write one folder per calculation, run them – here or on a
cluster – and read the results back. The settings come from an ASE ``Espresso`` calculator.
pw.x is installed separately (for example ``conda install -c conda-forge qe``), and the
pseudopotential file goes into the folder ``pseudo``.

**Step 1: settings and structure**

.. literalinclude:: scripts/al_qe.py
   :language: python
   :start-at: from ase.build import bulk
   :end-before: ec = ElasticConstants

**Step 2: write, run and read the run folders**

.. literalinclude:: scripts/al_qe.py
   :language: python
   :start-at: ec = ElasticConstants

``write_inputs`` creates ``runs/ground_state``, ``runs/C11_C12_I/strain_-0.010000`` and so on;
``run_jobs`` runs pw.x in every folder that has no usable output yet (call it again after a
failure, and only the missing runs are repeated); ``qe_status`` lists the runs that are done,
failed – with the reason – and pending. To run on a cluster instead, submit the folders
yourself and call ``read_results`` afterwards (:doc:`../guide/dft`).

**Screen output** (shortened):

.. literalinclude:: output/al_qe.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - 100.0
     - 68.5
     - 24.2

These settings (12³ k-points, 40 Ry) are a quick start, not converged for Al.
Script: :download:`al_qe.py <scripts/al_qe.py>`. VASP works the same way with ``code="vasp"``
and your INCAR, KPOINTS and POTCAR (:doc:`../guide/dft`).

Comparison
----------

.. list-table:: fcc aluminium, GPa
   :header-rows: 1
   :widths: 40 20 20 20

   * - Calculator
     - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - EAM (Mishin 1999)
     - 112.8
     - 60.5
     - 31.7
   * - GRACE-2L-OAM
     - 112.8
     - 60.5
     - 33.9
   * - MACE medium-mpa-0
     - 114.8
     - 64.9
     - 29.3
   * - GPAW, PBE, 32³ k-points
     - 118.8
     - 57.8
     - 33.9
   * - Quantum ESPRESSO, PBE, 12³ k-points
     - 100.0
     - 68.5
     - 24.2
   * - Experiment (room temperature)
     - 106.8
     - 60.7
     - 28.2

The differences are those of the calculators and their settings; the method is the same.
With any other ASE calculator, replace the calculator line – nothing else changes.
