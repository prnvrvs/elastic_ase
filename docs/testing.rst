Accuracy
========

How close are the elastic constants from elastic_ase to other methods and to experiment? The
method itself adds little error: the result is that of your calculator and its settings.

Same potential, other method
----------------------------

The other common way to compute elastic constants is the **stress–strain** method: apply
each strain component, compute the stress, and take :math:`C_{ij}` from the slope of the
stress. With the same potential (GRACE-2L-OAM), the same relaxed structures and the atoms
relaxed in every strained cell, both methods agree within about 1 GPa for every Laue group.
The small differences come from the two methods sampling the anharmonic (higher-order) terms
differently at the same strains.

Same potential, published constants
-----------------------------------

For the aluminium EAM potential of Mishin et al.
(`Phys. Rev. B 59, 3393 (1999) <https://doi.org/10.1103/PhysRevB.59.3393>`_, file
``Al99.eam.alloy`` from the `NIST Interatomic Potentials Repository
<https://www.ctcms.nist.gov/potentials/>`_):

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * -
     - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - Published
     - 114.0
     - 61.6
     - 31.6
   * - elastic_ase (strains ±0.25 %)
     - 113.5
     - 61.3
     - 31.6

The primitive cell, the conventional 4-atom cell and a rotated cell give identical constants.
With the default strains of ±1 %, :math:`C_{11}` is about 1 GPa lower: the anharmonic terms of
the potential start to contribute.

Same structure, any cell or orientation
---------------------------------------

The constants do not depend on how the crystal is given: for every Laue group, the default
route (rotation into the standard setting plus the named table), the full 21-constant fit and
the orientation-aware automatic mode on randomly rotated copies give the same stiffness tensor
to within a fraction of a GPa.

Calculated against measured
---------------------------

How close the constants are to experiment depends on the calculator. For GRACE-2L-OAM (trained
on PBE data) through elastic_ase, in GPa, against room-temperature measurements:

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

PBE and potentials trained on it typically underestimate the stiffness of strongly bonded
solids such as Si and MgO, and the measured values include thermal softening.

DFT
---

With a DFT code, the accuracy is set by its convergence – above all the k-point mesh for
metals. For fcc Al with PBE, the shear constants still change by up to 5 GPa between 28³ and 32³
k-points, while the total energy has been converged to 1 meV/atom since 20³. Converge the
elastic constants themselves (:doc:`guide/convergence`).

Results read from VASP OUTCARs agree with elastic_vasp, which uses the same deformation tables
(fcc Al with PBE: :math:`C_{11} = 106.6`, :math:`C_{12} = 63.4`, :math:`C_{44} = 28.9` GPa).

Limitations
-----------

* The constants are static (0 K, no zero-point motion); measured values at room temperature
  are usually a few percent lower.
* The symmetry analysis ignores magnetic order: for magnetic materials check that the
  magnetic structure does not lower the symmetry, or use ``symmetry="triclinic"``.
* The hardness models in the derived properties are empirical estimates.

Problems and questions: `issue tracker <https://github.com/prnvrvs/elastic_ase/issues>`_.
