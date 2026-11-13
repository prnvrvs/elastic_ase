Results and checks
==================

``calculate()``, ``read_results()`` and the ``vasp_output`` / ``qe_output`` shortcuts all
return an :class:`~elastic_ase.results.ElasticResults`, also kept as ``ec.results``. Its fields
can be read as attributes (``results.stiffness_matrix``) or as dictionary keys
(``results["stiffness_matrix"]``). Moduli are in GPa.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Key
     - Content
   * - ``stiffness_matrix``
     - 6×6 Voigt stiffness matrix in the frame chosen with ``frame=``
   * - ``stiffness_cartesian``, ``stiffness_standard``
     - the same in the calculation frame and in the standard crystallographic setting
   * - ``independent_constants``
     - the table's constants, e.g. ``{"C11": ..., "C12": ..., "C44": ...}``
   * - ``standard_constants``
     - the independent constants of the space group in the standard setting
   * - ``standard_rotation``
     - rotation from the calculation frame to the standard setting
   * - ``is_stable``, ``stability_criteria``
     - Born stability: all eigenvalues of C positive (the text names the smallest one)
   * - ``polycrystalline_averages``
     - see below
   * - ``derived_properties``
     - see below
   * - ``residual_stress``
     - ``pressure_GPa`` of the reference cell, ``max_along_deformation_GPa`` and
       ``per_deformation_GPa`` (see below)
   * - ``directional_properties``
     - with elastic-vis: ``E_min``, ``E_max``, ``E_anisotropy``, ``G_min``, ``G_max``,
       ``G_anisotropy``, ``nu_min``, ``nu_max``, ``is_auxetic``, ``LC_min``, ``LC_max``,
       ``is_negative_LC``; otherwise ``None``
   * - ``coefficients``, ``r2_scores``
     - fitted curvature :math:`a` of every deformation (GPa) and the :math:`R^2` of its fit
   * - ``energies``
     - every energy in eV, by run name
   * - ``fits``
     - strains, :math:`(E - E_0)/V_0` and the fitted polynomial of every deformation (used by
       ``plot_fits``)

Polycrystalline averages
------------------------

``K_voigt``, ``K_reuss``, ``K_VRH``, ``G_voigt``, ``G_reuss``, ``G_VRH`` (bulk and shear
moduli), ``E_VRH`` (Young's modulus), ``nu_VRH`` (Poisson's ratio), ``Pugh_ratio`` (:math:`K/G`;
above about 1.75 usually ductile), ``Cauchy_pressure`` (:math:`C_{12} - C_{44}`),
``Zener_anisotropy`` (:math:`2C_{44}/(C_{11} - C_{12})`, 1 for isotropic),
``universal_anisotropy`` (:math:`A^U = 5G_V/G_R + K_V/K_R - 6`, 0 for isotropic) and the
Chung–Buessem anisotropies ``chung_buessem_G``, ``chung_buessem_K``. Formulas in
:doc:`../theory`.

Derived properties
------------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Key
     - Meaning
   * - ``density_kg_m3``
     - mass density of the structure
   * - ``v_longitudinal_m_s``, ``v_transverse_m_s``, ``v_mean_m_s``
     - sound velocities of the isotropic polycrystal (Hill moduli)
   * - ``debye_temperature_K``
     - Debye temperature from the mean sound velocity (Anderson)
   * - ``kleinman_parameter``
     - cubic crystals only: :math:`\zeta = (C_{11} + 8C_{12})/(7C_{11} + 2C_{12})`
   * - ``hardness_chen_GPa``, ``hardness_tian_GPa``
     - empirical Vickers hardness models (Chen 2011, Tian 2012). Rough estimates fitted to
       hard covalent and ionic solids; for soft metals they are not meaningful.

For fcc Al with PBE (:math:`C_{11} = 106.6`, :math:`C_{12} = 63.4`, :math:`C_{44} = 28.9`
GPa) this gives :math:`v_l = 6422` m/s (experiment about 6420 m/s) and
:math:`\theta_D = 405` K (experiment 428 K).

The same function is available on its own:
``elastic_ase.derived_properties(averages, atoms, stiffness_standard, symmetry)``.

Residual stress
---------------

The method assumes a stress-free reference cell. If it is not – the cell was not relaxed, or
relaxed with other settings – the energy–strain curves pick up a linear term, and the
constants are biased by about the size of the stress (1.6 GPa of tension changes the
constants of EMT Al by about 0.8 GPa).

The linear coefficient :math:`b` of each fit equals :math:`\sigma_0 \cdot e` for that
deformation, so the fits themselves give the residual stress – no stress calculation is
needed, and it also works on OUTCARs. elastic_ase reports the residual pressure
:math:`P = -\operatorname{tr}\sigma_0/3` and warns above 0.5 GPa. The linear term is taken from
a cubic fit, because with a quadratic fit the anharmonic cubic term would show up as a fake
pressure at larger strains (0.3 GPa for Al at ±3 %).

.. code-block:: python

   results["residual_stress"]["pressure_GPa"]   # should be close to 0

Fit quality and plots
---------------------

.. code-block:: python

   ec.plot_fits("energy_strain.png")      # one panel per deformation, points + fitted parabola
   results["r2_scores"]                   # R^2 of every quadratic fit

Look at the curves first: scatter means unconverged settings (often the k-point mesh), a
single point off the curve a failed or mixed-up run, an asymmetric curve with a poor
:math:`R^2` strains that are too large.

Directional properties (elastic-vis)
------------------------------------

With elastic-vis installed:

.. code-block:: python

   ec.Young(theta, phi)                   # Young's modulus in a direction (GPa)
   ec.shear(theta, phi)                   # min and max shear modulus
   ec.poisson(theta, phi)                 # min and max Poisson's ratio
   ec.linear_compressibility(theta, phi)
   ec.visualize(properties=["Young", "Shear"], dimensions=["2D", "3D"],
                prefix="al", output_format="png")
   ec.to_elastic_vis()                    # elastic_vis.Elastic object

Saving and exporting
--------------------

.. code-block:: python

   ec.save("al.json")                       # structure, settings, every energy and the results
   ec = ElasticConstants.load("al.json")    # recomputes everything; nothing is run
   ec.export_matrix("C.txt")                # plain 6x6 matrix (numpy.loadtxt, elastic-vis)
   ec.compliance                            # S = C^-1 (1/GPa)

The JSON file keeps every energy, so the analysis can be redone later – for example with a
newer elastic_ase – without new calculations.
