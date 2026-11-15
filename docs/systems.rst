Crystal systems and deformations
================================

For each named table: the layout of the stiffness matrix in terms of the independent
constants, and every deformation – its run folder name, the deformation gradient
:math:`F = I + \varepsilon` applied to the cell (rows of lattice vectors :math:`L' = L F`) and
the combination of constants its fitted curvature :math:`a` measures. The tables are
generated from the code.

A term :math:`\frac{\delta^2}{k - \delta^2}` on the diagonal (:math:`k = 1` for opposite
normal strains, 4 for shears) keeps the cell volume constant (:math:`\det F = 1`). For a
stress-free cell it does not change :math:`a`; with a residual pressure it removes the
:math:`-P\,\Delta V` term (see :doc:`theory`).

.. include:: _generated/system_tables.rst
