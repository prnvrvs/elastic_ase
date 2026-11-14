Citing
======

If elastic_ase is useful for your work, please cite these papers (the same as for
`elastic_vasp <https://prnvrvs.github.io/elastic_vasp/>`_, which elastic_ase is based on):

1. P. Kumar and I. Adlakha, "Effect of Interstitial Hydrogen on Elastic Behavior of Metals:
   An Ab-Initio Study", *Journal of Engineering Materials and Technology* **145** (1), 2022.
   `doi:10.1115/1.4055097 <https://doi.org/10.1115/1.4055097>`_

2. P. Mishra, P. Kumar, L. Neelakantan and I. Adlakha, "First-principles prediction of
   electrochemical polarization and mechanical behavior in Mg based intermetallics",
   *Computational Materials Science* **214**, 111667, 2022.
   `doi:10.1016/j.commatsci.2022.111667 <https://doi.org/10.1016/j.commatsci.2022.111667>`_

BibTeX
------

.. code-block:: bibtex

   @article{Kumar_2022,
     author  = {Kumar, P. and Adlakha, I.},
     title   = {{Effect of Interstitial Hydrogen on Elastic
                 Behavior of Metals: An Ab-Initio Study}},
     journal = {Journal of Engineering Materials and Technology},
     volume  = {145},
     number  = {1},
     year    = {2022},
     doi     = {10.1115/1.4055097},
   }

   @article{Mishra_2022,
     author  = {Mishra, Pragyandipta and Kumar, Pranav and
                Neelakantan, Lakshman and Adlakha, Ilaksh},
     title   = {{First-principles prediction of electrochemical
                 polarization and mechanical behavior in Mg based
                 intermetallics}},
     journal = {Computational Materials Science},
     volume  = {214},
     pages   = {111667},
     year    = {2022},
     doi     = {10.1016/j.commatsci.2022.111667},
   }


Software elastic_ase relies on
------------------------------

Please also cite these as their authors request:

* `ASE <https://docs.ase-lib.org/>`_ – structures, calculators and file formats
  (`how to cite ASE <https://docs.ase-lib.org/faq.html#how-should-i-cite-ase>`_),
* `spglib <https://spglib.readthedocs.io/en/stable/>`_ – crystal symmetry (the papers to cite
  are listed at the top of its front page),
* `elastic-vis <https://github.com/prabhat7iitm/elastic-vis>`_ – if you use the directional
  properties or plots,
* the calculator or DFT code that produced the energies.

Formulas for the derived properties
-----------------------------------

If you report these quantities, the original references are (see :doc:`theory`):

* Voigt–Reuss–Hill averages: R. Hill, *Proc. Phys. Soc. A* **65**, 349 (1952).
* Universal anisotropy index: S. I. Ranganathan and M. Ostoja-Starzewski,
  *Phys. Rev. Lett.* **101**, 055504 (2008).
* Debye temperature from sound velocities: O. L. Anderson, *J. Phys. Chem. Solids* **24**,
  909 (1963).
* Hardness models: X.-Q. Chen et al., *Intermetallics* **19**, 1275 (2011); Y. Tian et al.,
  *Int. J. Refract. Met. Hard Mater.* **33**, 93 (2012).

Author and license
------------------

Pranav Kumar – `homepage <https://sites.google.com/view/prnvkmr4>`_. MIT license.
