Quick start
===========

A complete calculation with ASE's built-in EMT potential takes a few seconds:

.. code-block:: python

   from ase.build import bulk
   from ase.calculators.emt import EMT
   from elastic_ase import ElasticConstants

   # The structure must be relaxed with the same calculator: a pre-stressed cell biases the
   # constants. 3.9943 A is the EMT equilibrium lattice constant of Al.
   atoms = bulk("Al", "fcc", a=3.9943, cubic=True)
   atoms.calc = EMT()

   ec = ElasticConstants(atoms)          # symmetry detected: cubic
   results = ec.calculate()

   C = results["stiffness_matrix"]       # 6x6, GPa
   print(results["independent_constants"])   # {'C11': ..., 'C12': ..., 'C44': ...}
   print(results["polycrystalline_averages"]["K_VRH"], "GPa")
   print(results["is_stable"])

``calculate()`` prints a report: the fitted curvature of every deformation, the stiffness
matrix, Born stability, Voigt–Reuss–Hill moduli, derived properties (sound velocities, Debye
temperature, …) and, with elastic-vis installed, the directional extremes. Pass
``verbose=False`` to silence it.

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

.. code-block:: python

   ec.save("al_emt.json")                        # structure, settings, energies, results
   ec2 = ElasticConstants.load("al_emt.json")    # full analysis again, nothing is run
   ec.export_matrix("C.txt")                     # plain 6x6 matrix (numpy.loadtxt, elastic-vis)

A machine-learning potential
----------------------------

Only the calculator changes (not tested here – GRACE is not installed in our test setup):

.. code-block:: python

   from tensorpotential.calculator import grace_fm

   atoms = bulk("Co", "hcp", a=2.498, c=4.053)
   atoms.calc = grace_fm("GRACE-2L-OAM")
   # relax first (cell and atoms) with the same calculator, e.g. ase.filters.FrechetCellFilter
   results = ElasticConstants(atoms).calculate()

For multi-atom cells whose atoms move under strain, pass an ASE optimizer:
``ElasticConstants(atoms, optimizer=BFGS, fmax=1e-3)`` relaxes the atoms inside every strained
cell (:doc:`guide/calculators`).

Next steps
----------

* DFT with VASP or Quantum ESPRESSO: :doc:`guide/dft`.
* Your crystal is not cubic, or not in the standard orientation: :doc:`guide/orientation`.
* Everything in the results: :doc:`guide/results`.
