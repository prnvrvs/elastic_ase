Examples by crystal symmetry
============================

The crystal's Laue class decides how many independent elastic constants it has, and so which
deformations elastic_ase applies and how many calculations are needed. This page has one
example for each of the 11 Laue classes. Every example goes through the same four steps:

1. **Structure** – build or read it; ``check_symmetry`` shows which method will be used.
2. **Relax** – cell and atoms, with the calculator you use for the elastic constants.
3. **Elastic constants** – ``ElasticConstants(atoms).calculate()``; in multi-atom cells the
   atoms are relaxed inside every strained cell (``relaxation=True``, the default).
4. **Check and save** – the energy–strain curves and a JSON file with everything.

Steps 2 to 4 are the same for every crystal (shown here for silicon); each example below
shows its structure, the screen output and the result, with the complete script to download.

.. literalinclude:: scripts/si_diamond.py
   :language: python
   :start-after: # Step 2: relax

The calculator is the
machine-learning potential GRACE-2L-OAM (``pip install tensorpotential``), but any ASE
calculator works the same way. The numbers are those of this potential – they show what the
output looks like, not reference values for the materials. The Laue classes 6/m and -1 have no
simple common example, so a small model crystal shows their method.

.. list-table:: Overview
   :header-rows: 1
   :widths: 14 24 20 14 16

   * - Laue class
     - Example
     - Table
     - Constants
     - Calculations
   * - m-3m
     - Si
     - ``cubic``
     - 3
     - 13
   * - m-3
     - pyrite FeS₂
     - ``cubic``
     - 3
     - 13
   * - 6/mmm
     - Mg
     - ``hexagonal``
     - 5
     - 21
   * - 6/m
     - model crystal
     - ``hexagonal``
     - 5
     - 21
   * - -3m
     - corundum Al₂O₃
     - ``trigonal_2``
     - 7
     - 29
   * - -3
     - ilmenite FeTiO₃
     - ``trigonal_2``
     - 7
     - 29
   * - 4/mmm
     - rutile TiO₂
     - ``tetragonal_1``
     - 6
     - 25
   * - 4/m
     - scheelite CaWO₄
     - ``tetragonal_2``
     - 7
     - 29
   * - mmm
     - black phosphorus
     - ``orthogonal``
     - 9
     - 37
   * - 2/m
     - baddeleyite ZrO₂
     - ``monoclinic``
     - 13
     - 53
   * - -1
     - model crystal
     - ``triclinic``
     - 21
     - 85

"Calculations" counts the unstrained cell and every deformation at the default strains (±0.5 and ±1 %); the deformations are listed in :doc:`../systems`.

Cubic
-----

m-3m: silicon (diamond)
^^^^^^^^^^^^^^^^^^^^^^^

**Method:** ``cubic`` table – 3 constants (:math:`C_{11}, C_{12}, C_{44}`), 3 deformations, 13 calculations.

**Structure:**

.. literalinclude:: scripts/si_diamond.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/si_diamond.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - 142.8
     - 67.7
     - 76.2

:math:`K_{VRH} = 92.7` GPa, :math:`G_{VRH} = 57.4` GPa, mechanically stable. Script: :download:`si_diamond.py <scripts/si_diamond.py>`.

m-3: pyrite FeS₂
^^^^^^^^^^^^^^^^

**Method:** ``cubic`` table, as for m-3m: the lower point symmetry does not add elastic constants.

**Structure:**

.. literalinclude:: scripts/pyrite_fes2.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/pyrite_fes2.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{44}`
   * - 352.2
     - 17.9
     - 121.8

:math:`K_{VRH} = 129.3` GPa, :math:`G_{VRH} = 138.3` GPa, mechanically stable. Script: :download:`pyrite_fes2.py <scripts/pyrite_fes2.py>`.

Hexagonal
---------

6/mmm: magnesium (hcp)
^^^^^^^^^^^^^^^^^^^^^^

**Method:** ``hexagonal`` table – 5 constants (:math:`C_{11}, C_{12}, C_{13}, C_{33}, C_{44}`; :math:`C_{66} = (C_{11} - C_{12})/2`), 5 deformations, 21 calculations. The 6-fold axis is put along :math:`z` automatically.

**Structure:**

.. literalinclude:: scripts/mg_hcp.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/mg_hcp.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{33}`
     - :math:`C_{44}`
   * - 61.9
     - 27.2
     - 21.8
     - 65.8
     - 18.1

:math:`K_{VRH} = 36.8` GPa, :math:`G_{VRH} = 18.5` GPa, mechanically stable. Script: :download:`mg_hcp.py <scripts/mg_hcp.py>`.

6/m: model crystal (P6/m)
^^^^^^^^^^^^^^^^^^^^^^^^^

**Method:** ``hexagonal`` table, as for 6/mmm: both hexagonal Laue classes have the same 5 constants.

**Structure:**

.. literalinclude:: scripts/model_p6m.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Relax:** the cell is relaxed with the atoms kept in place – a model, not a real material.

.. literalinclude:: scripts/model_p6m.py
   :language: python
   :start-after: # Step 2: relax
   :end-before: # Step 3: elastic constants

**Elastic constants:** the atoms stay fixed in the strained cells as well.

.. literalinclude:: scripts/model_p6m.py
   :language: python
   :start-after: # Step 3: elastic constants
   :end-before: # Step 4: check and save

**Screen output** (shortened):

.. literalinclude:: output/model_p6m.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{33}`
     - :math:`C_{44}`
   * - 230.1
     - 115.1
     - -6.2
     - 623.4
     - 15.8

:math:`K_{VRH} = 138.1` GPa, :math:`G_{VRH} = 56.8` GPa, mechanically stable. Script: :download:`model_p6m.py <scripts/model_p6m.py>`.

Trigonal
--------

-3m: corundum Al₂O₃
^^^^^^^^^^^^^^^^^^^

**Method:** ``trigonal_2`` table – 7 constants (:math:`C_{11}, C_{12}, C_{13}, C_{33}, C_{44}, C_{14}, C_{15}`), 7 deformations, 29 calculations. For -3m the fitted :math:`C_{15}` comes out zero.

**Structure:**

.. literalinclude:: scripts/corundum_al2o3.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/corundum_al2o3.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{33}`
     - :math:`C_{44}`
     - :math:`C_{14}`
     - :math:`C_{15}`
   * - 439.0
     - 152.5
     - 109.3
     - 439.2
     - 125.5
     - 17.9
     - 0.0

:math:`K_{VRH} = 228.6` GPa, :math:`G_{VRH} = 139.9` GPa, mechanically stable. Script: :download:`corundum_al2o3.py <scripts/corundum_al2o3.py>`.

-3: ilmenite FeTiO₃
^^^^^^^^^^^^^^^^^^^

**Method:** ``trigonal_2`` table – 7 constants (:math:`C_{15}` in addition to those of -3m), 7 deformations, 29 calculations.

**Structure:**

.. literalinclude:: scripts/ilmenite_fetio3.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/ilmenite_fetio3.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{33}`
     - :math:`C_{44}`
     - :math:`C_{14}`
     - :math:`C_{15}`
   * - 258.4
     - 155.2
     - 104.5
     - 191.4
     - 5.7
     - 14.8
     - 6.2

:math:`K_{VRH} = 155.6` GPa, :math:`G_{VRH} = 18.5` GPa, mechanically stable. Script: :download:`ilmenite_fetio3.py <scripts/ilmenite_fetio3.py>`.

Tetragonal
----------

4/mmm: rutile TiO₂
^^^^^^^^^^^^^^^^^^

**Method:** ``tetragonal_1`` table – 6 constants (:math:`C_{11}, C_{12}, C_{13}, C_{33}, C_{44}, C_{66}`), 6 deformations, 25 calculations. Explained in detail in :doc:`step_by_step`.

**Structure:**

.. literalinclude:: scripts/rutile_tio2.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/rutile_tio2.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{33}`
     - :math:`C_{44}`
     - :math:`C_{66}`
   * - 262.0
     - 170.6
     - 144.1
     - 460.3
     - 103.5
     - 210.2

:math:`K_{VRH} = 207.1` GPa, :math:`G_{VRH} = 105.5` GPa, mechanically stable. Script: :download:`rutile_tio2.py <scripts/rutile_tio2.py>`.

4/m: scheelite CaWO₄
^^^^^^^^^^^^^^^^^^^^

**Method:** ``tetragonal_2`` table – 7 constants (:math:`C_{16}` in addition), 7 deformations, 29 calculations.

**Structure:**

.. literalinclude:: scripts/scheelite_cawo4.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/scheelite_cawo4.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{33}`
     - :math:`C_{44}`
     - :math:`C_{66}`
     - :math:`C_{16}`
   * - 123.5
     - 56.7
     - 43.1
     - 111.1
     - 22.7
     - 30.0
     - -17.6

:math:`K_{VRH} = 71.2` GPa, :math:`G_{VRH} = 26.8` GPa, mechanically stable. Script: :download:`scheelite_cawo4.py <scripts/scheelite_cawo4.py>`.

Orthorhombic
------------

mmm: black phosphorus
^^^^^^^^^^^^^^^^^^^^^

**Method:** ``orthogonal`` table – 9 constants (:math:`C_{11}, C_{22}, C_{33}, C_{12}, C_{13}, C_{23}, C_{44}, C_{55}, C_{66}`), 9 deformations, 37 calculations.

**Structure:**

.. literalinclude:: scripts/black_p.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/black_p.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{22}`
     - :math:`C_{33}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{23}`
     - :math:`C_{44}`
   * - 146.2
     - 7.3
     - 32.6
     - -3.2
     - 27.7
     - 0.5
     - 2.3
   * - :math:`C_{55}`
     - :math:`C_{66}`
     - 
     - 
     - 
     - 
     - 
   * - 35.1
     - 2.4
     - 
     - 
     - 
     - 
     - 

:math:`K_{VRH} = 16.1` GPa, :math:`G_{VRH} = 11.6` GPa, mechanically stable. Script: :download:`black_p.py <scripts/black_p.py>`.

Monoclinic
----------

2/m: baddeleyite ZrO₂
^^^^^^^^^^^^^^^^^^^^^

**Method:** ``monoclinic`` table – 13 constants (the 9 orthorhombic ones and :math:`C_{15}, C_{25}, C_{35}, C_{46}`), 13 deformations, 53 calculations. The unique axis :math:`b` is put along :math:`y` automatically.

**Structure:**

.. literalinclude:: scripts/baddeleyite_zro2.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Screen output** (shortened):

.. literalinclude:: output/baddeleyite_zro2.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{22}`
     - :math:`C_{33}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{23}`
     - :math:`C_{44}`
   * - 297.1
     - 347.7
     - 266.1
     - 152.4
     - 107.4
     - 142.0
     - 72.5
   * - :math:`C_{55}`
     - :math:`C_{66}`
     - :math:`C_{15}`
     - :math:`C_{25}`
     - :math:`C_{35}`
     - :math:`C_{46}`
     - 
   * - 81.0
     - 110.2
     - 38.6
     - -8.3
     - -0.7
     - -8.7
     - 

:math:`K_{VRH} = 186.7` GPa, :math:`G_{VRH} = 83.7` GPa, mechanically stable. Script: :download:`baddeleyite_zro2.py <scripts/baddeleyite_zro2.py>`.

Triclinic
---------

-1: model crystal (P-1)
^^^^^^^^^^^^^^^^^^^^^^^

**Method:** ``triclinic`` table – all 21 constants, 21 deformations, 85 calculations. It needs no symmetry, so it is also the fallback when the symmetry cannot be found, and a cross-check for any crystal (``symmetry="triclinic"``).

**Structure:**

.. literalinclude:: scripts/model_p1.py
   :language: python
   :start-after: # Step 1: structure
   :end-before: # Step 2: relax

**Relax:** the cell is relaxed with the atoms kept in place – a model, not a real material.

.. literalinclude:: scripts/model_p1.py
   :language: python
   :start-after: # Step 2: relax
   :end-before: # Step 3: elastic constants

**Elastic constants:** the atoms stay fixed in the strained cells as well.

.. literalinclude:: scripts/model_p1.py
   :language: python
   :start-after: # Step 3: elastic constants
   :end-before: # Step 4: check and save

**Screen output** (shortened):

.. literalinclude:: output/model_p1.txt
   :language: text

**Result:**

.. list-table:: Elastic constants (GPa)
   :class: ex-constants

   * - :math:`C_{11}`
     - :math:`C_{12}`
     - :math:`C_{13}`
     - :math:`C_{14}`
     - :math:`C_{15}`
     - :math:`C_{16}`
     - :math:`C_{22}`
   * - 193.6
     - 101.6
     - 76.8
     - 8.9
     - -0.4
     - 1.8
     - 150.9
   * - :math:`C_{23}`
     - :math:`C_{24}`
     - :math:`C_{25}`
     - :math:`C_{26}`
     - :math:`C_{33}`
     - :math:`C_{34}`
     - :math:`C_{35}`
   * - 103.1
     - 2.4
     - -12.3
     - 10.9
     - 182.9
     - 12.0
     - 5.7
   * - :math:`C_{36}`
     - :math:`C_{44}`
     - :math:`C_{45}`
     - :math:`C_{46}`
     - :math:`C_{55}`
     - :math:`C_{56}`
     - :math:`C_{66}`
   * - -1.3
     - 55.9
     - 0.7
     - -12.1
     - 23.2
     - 0.9
     - 44.1

:math:`K_{VRH} = 119.9` GPa, :math:`G_{VRH} = 36.2` GPa, mechanically stable. Script: :download:`model_p1.py <scripts/model_p1.py>`.

Structures in another orientation
---------------------------------

All examples start in the standard orientation. A structure in any other orientation gives the
same constants: by default elastic_ase rotates it rigidly into the standard setting (only the
orientation changes, not the cell or the atoms). See :doc:`../guide/orientation`.
