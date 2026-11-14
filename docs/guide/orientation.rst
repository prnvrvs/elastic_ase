Symmetry and orientation
========================

Which table is used
-------------------

With ``symmetry="auto"`` (the default) spglib determines the space group and the matching
named table is used:

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Space groups
     - Table
     - Independent constants
   * - 195–230 (cubic)
     - ``cubic``
     - 3
   * - 168–194 (hexagonal)
     - ``hexagonal``
     - 5
   * - 143–167 (trigonal)
     - ``trigonal_2``
     - 7 (fits :math:`C_{14}` and :math:`C_{15}`; covers 3m, 32 and -3m in any in-plane
       orientation and the -3 class)
   * - 89–142 (tetragonal, Laue 4/mmm)
     - ``tetragonal_1``
     - 6
   * - 75–88 (tetragonal, Laue 4/m)
     - ``tetragonal_2``
     - 7
   * - 16–74 (orthorhombic)
     - ``orthogonal``
     - 9
   * - 3–15 (monoclinic)
     - ``monoclinic``
     - 13 (unique axis b)
   * - 1–2 (triclinic)
     - ``triclinic``
     - 21

If spglib cannot determine the symmetry, the triclinic table (valid for every crystal) is used
with a warning. The symmetry analysis ignores magnetic order.

Check a structure without any calculation:

.. code-block:: bash

   elastic-ase POSCAR --check

.. code-block:: python

   from elastic_ase import check_symmetry
   report = check_symmetry("POSCAR")

The report gives the space group, the independent constants and which tables fit the crystal
in its orientation.

Orientation
-----------

The named tables assume the crystal axes along fixed Cartesian directions: cubic, tetragonal
and orthorhombic axes along :math:`x, y, z`; the 3-, 4- or 6-fold axis along :math:`z`; the
monoclinic unique axis :math:`b` along :math:`y`. Three ways to deal with a structure in
another orientation:

``reorient=True`` (default)
   The structure is rotated rigidly into the standard setting (the rotation spglib reports for
   the standardized cell). **Only the orientation changes**: your cell – primitive,
   conventional or supercell – the number and order of atoms and all per-atom data (magnetic
   moments, charges, tags, constraints) stay as given, so calculator settings such as the
   k-point mesh remain valid. A primitive rhombohedral cell stays primitive, with its 3-fold
   axis rotated onto :math:`z`.

``reorient=False`` with a named table
   The structure is used as given. elastic_ase checks *exactly* whether the table can measure
   this crystal in this orientation – whether every stiffness matrix the point group allows
   here can be represented by the table – and warns if not, because the constants would be
   wrong.

``reorient=False`` with ``symmetry="auto"``
   Orientation-aware mode: the symmetry-allowed stiffness matrices are derived from the
   crystal's rotations in its own orientation (a basis of the space of matrices with
   :math:`M C M^\mathsf{T} = C` for every point-group operation), and a minimal set of
   well-conditioned deformations is chosen for them. Works in any orientation.

Output frame
------------

``frame="cartesian"`` (default) returns C in the frame the calculation was done in (after any
reorientation); ``frame="standard"`` returns it in the standard crystallographic setting.
Both are always available as ``results["stiffness_cartesian"]`` and
``results["stiffness_standard"]``, and ``results["standard_constants"]`` gives the
independent constants in the standard setting. Stability, Zener anisotropy and Cauchy pressure
are always evaluated in the standard setting, where entries such as :math:`C_{12}` and
:math:`C_{44}` have their usual meaning; the Voigt–Reuss–Hill moduli do not depend on the
frame.

Common pitfalls
---------------

* **Trigonal "31m"-type groups** (P-31m, P-31c, P31m, P31c, P312, P3₁12, P3₂12): with
  :math:`a \parallel x` the 2-fold axes are perpendicular to :math:`a`, so the stiffness
  matrix has :math:`C_{15}` instead of :math:`C_{14}`. ``trigonal_2`` fits both, which is why
  automatic detection always uses it.
* **Monoclinic cells with unique axis c** (:math:`\gamma \ne 90°`): handled by
  ``reorient=True``; with ``reorient=False`` and ``symmetry="monoclinic"`` you get a warning.
* **Hexagonal cells with c along x**: handled by ``reorient=True``; warned otherwise.
