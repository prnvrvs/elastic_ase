Installation
============

Python 3.10 or newer is required; NumPy, ASE and spglib are installed automatically.

From PyPI:

.. code-block:: bash

   pip install elastic-ase

The current version from GitHub:

.. code-block:: bash

   pip install "git+https://github.com/prnvrvs/elastic_ase.git"

Optional: the directional plots and extremes (Young's modulus, shear modulus, Poisson's ratio
and linear compressibility surfaces) use
`elastic-vis <https://github.com/prabhat7iitm/elastic-vis>`_. It is not on PyPI, so install it
separately:

.. code-block:: bash

   pip install "git+https://github.com/prabhat7iitm/elastic-vis.git"

For development – an editable install with the tools to run the tests and build this
documentation:

.. code-block:: bash

   git clone https://github.com/prnvrvs/elastic_ase.git
   cd elastic_ase
   pip install -e ".[test,doc]"
   python -m unittest discover -s tests
   make -C docs html

Calculators
-----------

elastic_ase needs only an ASE calculator that returns the energy; install the one you use
separately. Tested here: EMT and EAM (built into ASE), GRACE, MACE, GPAW and VASP.
Other machine-learning potentials (CHGNet, …) and ASE calculators should work the same way but
have not been tested yet (see :doc:`testing`).

Check the installation:

.. code-block:: bash

   elastic-ase --help
   python -c "import elastic_ase; print(elastic_ase.__version__)"
