Installation
============

Python 3.10 or newer is required; NumPy, ASE and spglib are installed automatically.

From PyPI:

.. code-block:: bash

   pip install elastic-ase

The current version from GitHub:

.. code-block:: bash

   pip install "git+https://github.com/prnvrvs/elastic_ase.git"

Calculators
-----------

elastic_ase needs an ASE calculator that returns the energy (and the forces, to relax the
atoms in multi-atom cells); install the one you use separately. EMT and EAM are built into ASE;
machine-learning potentials (GRACE, MACE, …) and GPAW are installed on their own
(:doc:`examples/calculators`), as are VASP and QE (:doc:`guide/dft`).

Check the installation
----------------------

.. code-block:: bash

   elastic-ase --help
   python -c "import elastic_ase; print(elastic_ase.__version__)"

Optional: directional properties (elastic-vis)
----------------------------------------------

The directional plots and extremes – Young's modulus, shear modulus, Poisson's ratio and
linear compressibility as a function of direction – use
`elastic-vis <https://github.com/prabhat7iitm/elastic-vis>`_. It is not on PyPI, so it cannot
be installed together with elastic_ase; install it separately:

.. code-block:: bash

   pip install "git+https://github.com/prabhat7iitm/elastic-vis.git"

With elastic-vis installed, the report printed by ``calculate()`` ends with one more block
(here for EMT aluminium):

.. code-block:: text

   ----------------------------------------------------------------------
    Directional Elastic Extremes & Anisotropy (elastic-vis):
     Young's Modulus E    : min =    28.21 GPa, max =    83.21 GPa  (Anisotropy: 2.95)
     Shear Modulus G      : min =    10.21 GPa, max =    36.18 GPa  (Anisotropy: 3.54)
     Poisson's Ratio nu   : min =   -0.227,     max =    0.756      (Negative nu: YES (auxetic))
     Linear Compress. LC  : min =     8.42,     max =     8.42 TPa⁻¹ (Negative LC: No)

Without elastic-vis everything else works – the stiffness matrix, independent constants,
stability, Voigt–Reuss–Hill averages and derived properties. Only these need it:

* ``results.directional_properties`` (``None`` without elastic-vis);
* ``ec.Young``, ``ec.shear``, ``ec.poisson``, ``ec.linear_compressibility``,
  ``ec.visualize`` and ``ec.to_elastic_vis`` (:doc:`guide/results`);
* the ``--plot`` options of the :doc:`command line <guide/cli>`.

Development
-----------

An editable install with the tools to run the tests and build this
documentation:

.. code-block:: bash

   git clone https://github.com/prnvrvs/elastic_ase.git
   cd elastic_ase
   pip install -e ".[test,doc]"
   python -m unittest discover -s tests
   make -C docs html
