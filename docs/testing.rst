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

**Every Laue class with GRACE-2L-OAM** (structure relaxed with its symmetry fixed; clamped
ions under strain; the three routes compared through the frame-independent eigenvalues of C)

.. list-table::
   :header-rows: 1
   :widths: 10 28 16 14 16 16

   * - Laue
     - Structure (space group)
     - Table
     - Constants
     - default vs triclinic
     - auto, rotated
   * - m-3m
     - NaCl (225)
     - ``cubic``
     - 3
     - 0.02 GPa
     - 0.04 GPa
   * - m-3
     - pyrite FeS₂ (205)
     - ``cubic``
     - 3
     - 0.87 GPa
     - 0.70 GPa
   * - 6/mmm
     - Mg (194)
     - ``hexagonal``
     - 5
     - 0.05 GPa
     - 0.06 GPa
   * - 6/m
     - model crystal (175)
     - ``hexagonal``
     - 5
     - 0.85 GPa
     - 0.48 GPa
   * - -3m
     - corundum Al₂O₃ (167)
     - ``trigonal_2``
     - 6
     - 0.12 GPa
     - 0.23 GPa
   * - -3
     - ilmenite FeTiO₃ (148)
     - ``trigonal_2``
     - 7
     - 0.21 GPa
     - 0.39 GPa
   * - 4/mmm
     - rutile TiO₂ (136)
     - ``tetragonal_1``
     - 6
     - 0.42 GPa
     - 0.40 GPa
   * - 4/m
     - scheelite CaWO₄ (88)
     - ``tetragonal_2``
     - 7
     - 0.08 GPa
     - 0.11 GPa
   * - mmm
     - model crystal (47)
     - ``orthogonal``
     - 9
     - 0.18 GPa
     - 0.14 GPa
   * - 2/m
     - model crystal (10)
     - ``monoclinic``
     - 13
     - 0.05 GPa
     - 0.14 GPa
   * - -1
     - model crystal (2)
     - ``triclinic``
     - 21
     - 0.25 GPa
     - 0.25 GPa

The model crystals (cell relaxed, atoms fixed) only check the method; their constants are not
physical. The remaining differences come from the different strains of each route probing the
anharmonic terms. This run also found that the Laue class reported for space groups 168–176
(6/m) and 195–206 (m-3) was wrong (fixed; tables and constants were not affected).

**Energy–strain (elastic_ase) against stress–strain (the Materials Project method)**, both with
GRACE-2L-OAM on the same relaxed structure, atoms relaxed in every strained cell (stress–strain:
each strain component at ±0.5 and ±1 %, :math:`C_{ij}` from the slope of the stress). Largest
difference over all entries of C:

.. list-table::
   :header-rows: 1
   :widths: 34 14 12 16 24

   * - Material
     - Laue
     - Atoms
     - Constants
     - max difference
   * - Al
     - m-3m
     - 1
     - 3
     - 0.05 GPa
   * - NaCl
     - m-3m
     - 2
     - 3
     - 0.03 GPa
   * - MgO
     - m-3m
     - 2
     - 3
     - 0.13 GPa
   * - Cu
     - m-3m
     - 1
     - 3
     - 0.10 GPa
   * - Si
     - m-3m
     - 2
     - 3
     - 0.41 GPa
   * - pyrite FeS₂
     - m-3
     - 12
     - 3
     - 0.95 GPa
   * - Mg
     - 6/mmm
     - 2
     - 5
     - 0.05 GPa
   * - Ti
     - 6/mmm
     - 2
     - 5
     - 0.37 GPa
   * - ZnO (wurtzite)
     - 6/mmm
     - 4
     - 5
     - 0.54 GPa
   * - rutile TiO₂
     - 4/mmm
     - 6
     - 6
     - 0.98 GPa
   * - scheelite CaWO₄
     - 4/m
     - 12
     - 7
     - 0.15 GPa
   * - corundum Al₂O₃
     - -3m
     - 10
     - 6
     - 0.19 GPa
   * - ilmenite FeTiO₃
     - -3
     - 10
     - 7
     - 1.21 GPa
   * - black P
     - mmm
     - 4
     - 9
     - 0.04 GPa
   * - ZrO₂ (baddeleyite)
     - 2/m
     - 12
     - 13
     - 2.07 GPa

The remaining differences (≤ 1 GPa except 2 GPa for one entry of ZrO₂) come from the two
methods probing anharmonic terms differently at the same strains.

**Against experiment** (same GRACE calculations; GPa; experimental values at room temperature
from the references given). The differences are those of the potential, which was trained
on PBE data – PBE typically underestimates the stiffness of Si and MgO – not of elastic_ase:

.. list-table::
   :header-rows: 1
   :widths: 14 43 43

   * - Material
     - GRACE through elastic_ase
     - Experiment
   * - Al
     - :math:`C_{11}` 112.8, :math:`C_{12}` 60.5, :math:`C_{44}` 33.9
     - 106.8, 60.7, 28.2 (Kamm & Alers, J. Appl. Phys. 35, 327, 1964)
   * - Cu
     - 175.1, 125.3, 76.4
     - 170.0, 122.5, 75.8 (Chang & Himmel, J. Appl. Phys., 1966)
   * - Si
     - 142.8, 67.7, 76.2
     - 165.7, 63.9, 79.6 (Mason 1958, as compiled in J. Synchrotron Rad. 2014)
   * - MgO
     - 258.1, 94.4, 128.3
     - 297.8, 95–97, 155.8 (inelastic X-ray scattering)
   * - Mg
     - :math:`C_{11}, C_{12}, C_{13}, C_{33}, C_{44}`: 61.9, 27.2, 21.8, 65.8, 18.1
     - 59.4, 25.6, 21.4, 61.6, 16.4 (Slutsky & Garland, Phys. Rev. 107, 972, 1957)
   * - Ti
     - 137.1, 92.8, 79.6, 181.3, 42.5
     - 162.4, 92.0, 69.0, 180.7, 46.7 (Fisher & Renken, Phys. Rev. 135, A482, 1964)
   * - rutile
     - :math:`C_{11}, C_{12}, C_{13}, C_{33}, C_{44}, C_{66}`: 262.0, 170.6, 144.1, 460.3, 103.5, 210.2
     - 268, 175, 147, 484, 123.8, 190.2 (Carnes, as quoted in arXiv:2607.08860)

A comparison with the Materials Project's own (PBE, stress–strain) values is not included: it
needs their API key.

**Other ASE calculators** (1-atom fcc Al, cell relaxed with the same calculator first)

* GRACE-2L-OAM (machine-learning potential): :math:`a_0 = 4.043` Å, :math:`C_{11} = 112.8`,
  :math:`C_{12} = 60.5`, :math:`C_{44} = 33.9` GPa, stable, residual pressure from the fits
  0.000 GPa; the command line (``-c grace``) gives the same numbers.
* GPAW (PBE, plane waves, also with MPI): the route works – relaxation, 13 runs, analysis, the
  primitive and the conventional cell – and the bulk modulus (76–80 GPa) agrees with VASP for
  every setting. The shear constants of aluminium, however, depend strongly on the k-mesh
  (Fermi–Dirac and Methfessel–Paxton smearing, 12³–30³): the cutoff (400, 600, 800 eV), the SCF
  tolerance and the pseudopotential (GPAW's PAW and its built-in HGH) change nothing, but
  :math:`C_{11} - C_{12}` from the (1, −1, 0) deformation alone is −2.0 GPa at 20³ and 64.7 GPa
  at 30³ – its whole energy curve spans only 0.66 meV at ±1 % strain. VASP with the tetrahedron
  method is not fully converged at these meshes either (:math:`C_{44}` between 24 and 36 GPa
  from 16³ to 32³, while its total energy is within 1 meV/atom from 20³ on). The
  residual-pressure check cannot detect this; see :doc:`guide/convergence`.
* MACE (``medium-mpa-0``, double precision): :math:`a_0 = 4.039` Å, :math:`C_{11} = 114.8`,
  :math:`C_{12} = 64.9`, :math:`C_{44} = 29.3` GPa, stable; the command line
  (``-c mace.calculators:mace_mp``) gives the same within 0.2 GPa (single precision). (On macOS,
  a conda NumPy together with a pip PyTorch can load two OpenMP runtimes and abort on import –
  an installation problem outside elastic_ase.)

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
* Machine-learning potentials other than GRACE and MACE (CHGNet, …): they should work like any
  ASE calculator, but have not been run here.
* GPAW with converged settings; other DFT codes (CASTEP, ABINIT, …) through
  ``write_inputs`` / ``read_results``.
* Magnetic materials (the symmetry analysis ignores magnetic order), and ``optimizer=`` with
  DFT calculators.
* The hardness models are empirical; they were only evaluated, not validated.
