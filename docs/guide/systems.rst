Crystal systems and deformations
================================

For each named table: the layout of the stiffness matrix in terms of the independent
constants, and every deformation – its run folder name, the deformation gradient
:math:`F = I + \varepsilon` applied to the cell (rows of lattice vectors :math:`L' = L F`) and
the combination of constants its fitted curvature :math:`a` measures. These tables are
generated from ``elastic_ase.systems`` when the documentation is built, so they always match
the code.

A term :math:`\frac{\delta^2}{k - \delta^2}` on the diagonal keeps the cell volume constant
(:math:`\det F = 1`); it is of second order and does not change :math:`a` (see
:doc:`../theory`).

.. include:: ../_generated/system_tables.rst
