What is tested
==============

This page lists what has been checked, and how – and what has not. Please report problems on
the `issue tracker <https://github.com/prnvrvs/elastic_ase/issues>`_.

Checked
-------

**Deformation tables and solver**

* Every named table recovers a known stiffness matrix from a synthetic quadratic energy
  (:math:`E = \tfrac12 V_0\,\boldsymbol\varepsilon^{\mathsf T} C \boldsymbol\varepsilon` with
  the Green–Lagrange strain) to within 0.01 GPa.
* The tables (constants, matrix layout, deformations) are identical to elastic_vasp's.
* Voigt–Reuss–Hill averages agree with elastic_vasp's.

**Symmetry, orientation, consistency (EMT)**

* Relaxed EMT alloy structures of every crystal family – cubic (Pm-3m), hexagonal
  (P6₃/mmc), tetragonal (P4/mmm, P4/m), trigonal (R-3m in the hexagonal and in the primitive
  rhombohedral cell, P-3), orthorhombic, monoclinic (unique axis b and c) and triclinic – give
  the same constants (to 0.2 GPa) from the default route (rotation into the standard setting
  plus the named table), from the full 21-constant triclinic fit and from the automatic mode
  on randomly rotated copies.
* A misoriented named table is detected (orthorhombic and monoclinic examples), the
  stability check is independent of the orientation, and the rotation keeps cell, atoms,
  magnetic moments, tags and constraints.

**VASP**

* Reading real VASP 6 OUTCARs: the 19 runs of fcc Al from the
  `elastic_vasp tutorial <https://prnvrvs.github.io/elastic_vasp/>`_, arranged as elastic_ase
  run folders, give :math:`C_{11} = 106.6`, :math:`C_{12} = 63.4`, :math:`C_{44} = 28.9` GPa,
  as elastic_vasp; a swapped OUTCAR is detected by the cell check.
* The electronic-convergence check gives the expected answer on ASE's sample OUTCARs.
* Writing run folders with your own INCAR, KPOINTS and POTCAR; the POTCAR element-order check
  (with a POTCAR header written for the test).
* The residual pressure from the fits of the real Al runs is 0.01 GPa (the cell is at the
  Birch–Murnaghan minimum).

**Quantum ESPRESSO**

* ``espresso.pwi`` written from an ASE ``Espresso`` calculator contains its settings
  (cutoff, k-points, pseudopotentials).
* Reading energies and the SCF-convergence check on the sample pw.x outputs shipped with ASE.

**Run management**

* Saved energies, failed runs, resuming only missing runs, ignoring energies saved for another
  structure: tested with a stand-in file-based calculator (EMT with a ``directory`` that can be
  made to fail or report an unconverged SCF).
* ``run_jobs``: tested with a stand-in program (a Python script that reads the POSCAR, computes
  the EMT energy and can be made to crash).
* Output folders with missing, unreadable or mixed-up files are reported.

**Derived properties**

* For fcc Al with the PBE constants above: :math:`v_l = 6422` m/s (experiment about
  6420 m/s), :math:`\theta_D = 405` K (experiment 428 K).
* Residual pressure from the fits equals EMT's stress (−1.561 GPa for Al at 4.05 Å).

The unit tests (``python -m unittest discover -s tests``, 47 tests) cover most of the above
with EMT and small files; the real OUTCARs are not part of the repository.

Not tested yet
--------------

* A real pw.x (Quantum ESPRESSO) run through elastic_ase, from input to result.
* ``run_jobs`` or ``calculate(workdir=...)`` with a real DFT program, and ASE's ``Vasp``
  calculator writing complete inputs (needs a POTCAR library).
* Machine-learning potentials (GRACE, MACE, CHGNet, …): they should work like any ASE
  calculator, but have not been run here.
* Other DFT codes (CASTEP, ABINIT, GPAW, …) through ``write_inputs`` / ``read_results``.
* Magnetic materials (the symmetry analysis ignores magnetic order), and ``optimizer=`` with
  DFT calculators.
* The hardness models are empirical; they were only evaluated, not validated.
