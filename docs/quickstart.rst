Quick start
===========

A complete calculation with ASE's built-in EMT potential takes a few seconds. The structure
must be relaxed with the same calculator – a pre-stressed cell biases the constants – so the
example uses 3.9943 Å, the EMT equilibrium lattice constant of Al:

.. code-block:: python

   from ase.build import bulk
   from ase.calculators.emt import EMT
   from elastic_ase import ElasticConstants

   atoms = bulk("Al", "fcc", a=3.9943, cubic=True)
   atoms.calc = EMT()

   ec = ElasticConstants(atoms)
   results = ec.calculate()

   C = results["stiffness_matrix"]       # 6x6, GPa
   print(results["independent_constants"])
   print(results["polycrystalline_averages"]["K_VRH"], "GPa")
   print(results["is_stable"])

Output on the screen
--------------------

``calculate()`` prints a report while it works:

.. code-block:: text

   [ElasticASE] Auto-detected crystal symmetry: 'cubic'
   ======================================================================
    ElasticASE: Elastic Constants Calculation (CUBIC)
   ======================================================================
   Initial Cell Volume V0 : 63.7268 Å³ (15.9317 Å³/atom)
   Strain values          : [-0.01, -0.005, 0.0, 0.005, 0.01]
   Deformation modes      : 3
   Calculations           : 13 (including the unstrained cell)
   Atom relaxation        : not needed (1 atom)
   ----------------------------------------------------------------------
     Mode C11_C12_I          -> Coeff =      20.43 GPa  (R² = 1.000000)
     Mode C11_C12_II         -> Coeff =     178.23 GPa  (R² = 0.999954)
     Mode C44                -> Coeff =      54.27 GPa  (R² = 0.999983)

   ======================================================================
    6x6 Stiffness Matrix C_ij (GPa) [Cartesian Frame]:
   ======================================================================
          C1     C2     C3     C4     C5     C6
   C1  53.23  32.80  32.80   0.00   0.00   0.00
   C2  32.80  53.23  32.80   0.00   0.00   0.00
   C3  32.80  32.80  53.23   0.00   0.00   0.00
   C4   0.00   0.00   0.00  36.18   0.00   0.00
   C5   0.00   0.00   0.00   0.00  36.18   0.00
   C6   0.00   0.00   0.00   0.00   0.00  36.18

   ----------------------------------------------------------------------
    Born Mechanical Stability:
     Criteria : C11 - C12 > 0 and C11 + 2*C12 > 0 and C44 > 0 (smallest eigenvalue 20.43 GPa)
     Status   : MECHANICALLY STABLE

   ----------------------------------------------------------------------
    Polycrystalline Moduli (Voigt-Reuss-Hill):
     Bulk Modulus K_VRH   :    39.61 GPa  (Voigt: 39.61, Reuss: 39.61)
     Shear Modulus G_VRH  :    21.87 GPa  (Voigt: 25.79, Reuss: 17.94)
     Young's Modulus E    :    55.40 GPa
     Poisson's Ratio nu   :    0.267
     Pugh Ratio (K/G)     :     1.81  (> 1.75 usually ductile)
     Cauchy Pressure      :    -3.38 GPa
     Zener Anisotropy A   :     3.54
     Universal Anisotropy : A^U =  2.190  (0 for isotropic)

   ----------------------------------------------------------------------
    Derived Properties:
     Density              :   2812.2 kg/m³
     Sound velocity       : longitudinal 4945 m/s, transverse 2788 m/s, mean 3102 m/s
     Debye temperature    :    367.0 K  (Anderson)
     Kleinman parameter   :    0.720
     Vickers hardness     : 3.07 GPa (Chen), 4.16 GPa (Tian)  [empirical, rough]
     Residual pressure    :    -0.00 GPa  (from the fits; should be ~0)
   ======================================================================

followed by the three ``print`` lines:

.. code-block:: text

   {'C11': 53.22540646083159, 'C12': 32.79886407889344, 'C44': 36.182294384609484}
   39.607711539539494 GPa
   True

How to read the report, from top to bottom:

* **Setup** – the symmetry found, the strains, and the number of calculations.
* **Mode … Coeff** – the curvature of each energy–strain curve and the :math:`R^2` of its
  fit; an :math:`R^2` well below 0.999 points to noisy energies.
* **Stiffness matrix** – the 6×6 :math:`C_{ij}` in GPa.
* **Born stability** – whether the crystal is mechanically stable.
* **Polycrystalline moduli** – bulk, shear and Young's modulus, Poisson's ratio and
  anisotropy indices of the isotropic polycrystal.
* **Derived properties** – density, sound velocities, Debye temperature, hardness estimates,
  and the **residual pressure**, which should be close to zero.
* With elastic-vis installed, one more block gives the directional extremes
  (:doc:`installation`).

``verbose=False`` switches the report off; the results are the same.

What happened
-------------

1. spglib found space group Fm-3m, so the **cubic** table is used: three deformations,
   each applied at the strains −1 %, −0.5 %, +0.5 % and +1 % (the default; change it with
   ``strains=[...]``). With the unstrained cell that makes 13 energies.
2. For every deformation the energy per volume is fitted with
   :math:`(E - E_0)/V_0 = a\,\delta^2 + b\,\delta + c`; the curvature :math:`a` is a known
   combination of elastic constants (e.g. :math:`a = C_{11} - C_{12}`).
3. The three curvatures are solved for :math:`C_{11}, C_{12}, C_{44}` (see :doc:`theory`).

Look at the curves
------------------

.. code-block:: python

   ec.plot_fits("energy_strain.png")

Smooth parabolas through all points are the first check. Scatter or kinks mean unconverged
settings or a failed calculation.

Save the results
----------------

``save`` keeps the structure, settings, every energy and the results; ``load`` redoes the
whole analysis from them without running anything; ``export_matrix`` writes the plain 6×6
matrix (for ``numpy.loadtxt`` or elastic-vis):

.. code-block:: python

   ec.save("al_emt.json")
   ec2 = ElasticConstants.load("al_emt.json")
   ec.export_matrix("C.txt")

Next steps
----------

* A complete calculation for a real crystal, step by step: :doc:`examples/step_by_step`.
* Your crystal symmetry: :doc:`examples/laue`. Your calculator: :doc:`examples/calculators`.
* DFT with VASP or Quantum ESPRESSO: :doc:`guide/dft`.
* Your crystal is not cubic, or not in the standard orientation: :doc:`guide/orientation`.
* Everything in the results: :doc:`guide/results`.
