Installation
============

Python 3.10 or newer is required; NumPy, ASE and spglib are installed automatically.

.. code-block:: bash

   pip install elastic-ase                                        # from PyPI (once the first release is published)
   pip install "git+https://github.com/prnvrvs/elastic_ase.git"   # current version from GitHub

Optional: the directional plots and extremes (Young's modulus, shear modulus, Poisson's ratio
and linear compressibility surfaces) use
`elastic-vis <https://github.com/prabhat7iitm/elastic-vis>`_. It is not on PyPI, so install it
separately:

.. code-block:: bash

   pip install "git+https://github.com/prabhat7iitm/elastic-vis.git"

For development (editable install with tests and documentation tools):

.. code-block:: bash

   git clone https://github.com/prnvrvs/elastic_ase.git
   cd elastic_ase
   pip install -e ".[test,doc]"
   python -m unittest discover -s tests     # run the tests
   make -C docs html                        # build this documentation

Calculators
-----------

elastic_ase needs only an ASE calculator that returns the energy; install the one you use
separately. Tested here: EMT (built into ASE) and VASP output files. Machine-learning
potentials (GRACE, MACE, CHGNet, …) and other ASE calculators should work the same way but
have not been tested yet (see :doc:`testing`).

Check the installation:

.. code-block:: bash

   elastic-ase --help
   python -c "import elastic_ase; print(elastic_ase.__version__)"
