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

**Running VASP through elastic_ase** (VASP 6.4.2 on macOS, 1-atom fcc Al, ENCUT 400 eV and a
7×7×7 k-mesh – deliberately small and *not converged* for Al, so this checks that every route
works, not the values):

* ``write_inputs`` → ``run_jobs`` → ``read_results``: 13 real runs, cubic symmetry detected.
* A randomly rotated copy of the cell with the default rotation into the standard setting
  gives the same constants (to 0.001 GPa).
* A crashed run (POTCAR removed) and an unconverged run (``NELM = 3``) are reported as failed
  with their reasons; ``read_results`` refuses to fit; ``write_inputs(runs=missing)`` and a
  second ``run_jobs`` repeat only those two runs, and the result equals the uninterrupted one.
* ASE's ``Vasp`` calculator: ``write_inputs`` writes complete inputs (INCAR, KPOINTS, POSCAR,
  POTCAR) from it, and ``calculate(workdir=...)`` runs VASP in one folder per run and, after
  being stopped part-way, resumes with only the missing runs – same constants as the run
  folders.
* The orientation-aware automatic mode on the rotated cell: VASP itself stopped in 2 of 13
  strained cells (its symmetry analysis found the direct and reciprocal lattices inconsistent,
  "I REFUSE TO CONTINUE WITH THIS SICK JOB"); elastic_ase reported both as failed. Set
  ``ISYM = 0`` (or change ``SYMPREC``) in the INCAR for rotated or low-symmetry cells; the
  default rotation into the standard setting avoids it.

**Against published values: Al EAM potential of Mishin et al.**
(`Phys. Rev. B 59, 3393 (1999) <https://doi.org/10.1103/PhysRevB.59.3393>`_; file
``Al99.eam.alloy`` from the `NIST Interatomic Potentials Repository
<https://www.ctcms.nist.gov/potentials/>`_, ASE's ``EAM`` calculator)

* The energy minimum is at 4.0500 Å, as published; residual pressure from the fits 0.000 GPa.
* With strains of ±0.25 %: :math:`C_{11} = 113.5`, :math:`C_{12} = 61.3`,
  :math:`C_{44} = 31.6` GPa – published 114.0, 61.6, 31.6 GPa.
* The primitive cell, the conventional 4-atom cell and a rotated cell (rotated back by
  default) give identical constants; the full triclinic fit and the automatic mode on the
  rotated cell agree to 1 GPa. With the default strains (±1 %) the anharmonic terms of this
  potential lower :math:`C_{11}` by about 1 GPa, slightly differently for each set of
  deformations.
* Debye temperature 431 K (experiment 428 K).

**Other ASE calculators** (1-atom fcc Al, cell relaxed with the same calculator first)

* GRACE-2L-OAM (machine-learning potential): :math:`a_0 = 4.043` Å, :math:`C_{11} = 112.8`,
  :math:`C_{12} = 60.5`, :math:`C_{44} = 33.9` GPa, stable, residual pressure from the fits
  0.000 GPa; the command line (``-c grace``) gives the same numbers.
* GPAW (PBE, plane waves 400 eV, Fermi–Dirac 0.1 eV): the route works (relaxation, 13 runs,
  analysis), but the values depend strongly on the k-mesh. With 8×8×8 the shear energies are
  noise (:math:`C_{44} < 0`); the curvature of the :math:`C_{44}` deformation gives
  −4, 52, 45 and 33 GPa at 8³, 12³, 16³ and 20³; at 20³ the bulk modulus (76.5 GPa) agrees with
  VASP (77.8 GPa) but :math:`C_{11} - C_{12}` (15 GPa) does not (VASP 43 GPa). Not converged –
  a reminder that metals need dense k-meshes, which the residual-pressure check cannot detect.
* MACE could not be tested: the local MACE installation fails on import (two OpenMP runtimes).

**Directional properties (elastic-vis)**

* All eight plots (Young's modulus, shear modulus, Poisson's ratio, linear compressibility;
  2D and 3D) are written; :math:`E_{[001]} = 1/S_{11}` and the shear modulus along [001]
  equals :math:`C_{44}` for cubic EMT Al.

**Quantum ESPRESSO**

* ``espresso.pwi`` written from an ASE ``Espresso`` calculator contains its settings
  (cutoff, k-points, pseudopotentials).
* Reading energies and the SCF-convergence check on the sample pw.x outputs shipped with ASE.

**Run management**

* Saved energies, failed runs, resuming only missing runs, ignoring energies saved for another
  structure: tested with real VASP (above) and, in the unit tests, with a stand-in file-based
  calculator (EMT with a ``directory`` that can be made to fail or report an unconverged SCF).
* ``run_jobs``: tested with real VASP (above) and, in the unit tests, with a stand-in program.
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
* VASP with converged settings through ``run_jobs`` (the VASP values above come from earlier
  runs read back with ``read_results``), larger cells, and VASP in parallel (MPI).
* Machine-learning potentials other than GRACE (MACE, CHGNet, …): they should work like any
  ASE calculator, but have not been run here.
* GPAW with converged settings; other DFT codes (CASTEP, ABINIT, …) through
  ``write_inputs`` / ``read_results``.
* Magnetic materials (the symmetry analysis ignores magnetic order), and ``optimizer=`` with
  DFT calculators.
* The hardness models are empirical; they were only evaluated, not validated.
