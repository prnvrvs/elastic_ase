Theory
======

Energy–strain relation
----------------------

For small homogeneous strains of a stress-free crystal the energy per reference volume is

.. math::

   \frac{E(\boldsymbol\varepsilon) - E_0}{V_0} = \frac12 \sum_{i,j=1}^{6} C_{ij}\,\varepsilon_i\,\varepsilon_j
   + \mathcal{O}(\varepsilon^3),

with the strain in Voigt notation, :math:`\varepsilon_1 = \varepsilon_{xx}`, …,
:math:`\varepsilon_3 = \varepsilon_{zz}` and engineering shear strains
:math:`\varepsilon_4 = 2\varepsilon_{yz}`, :math:`\varepsilon_5 = 2\varepsilon_{xz}`,
:math:`\varepsilon_6 = 2\varepsilon_{xy}`.

Each deformation uses a fixed strain pattern :math:`\mathbf e` scaled by the strain values
:math:`\delta`, :math:`\boldsymbol\varepsilon = \delta\,\mathbf e`, so

.. math::

   \frac{E(\delta) - E_0}{V_0} = a\, \delta^2 + b\,\delta + c + \mathcal{O}(\delta^3), \qquad
   a = \tfrac12\, \mathbf e^{\mathsf T} C\, \mathbf e .

Each curve is fitted with a quadratic; only the curvature :math:`a` is used for the constants
(converted from eV/Å³ to GPa with 1 eV/Å³ = 160.2177 GPa). :math:`b` absorbs a residual stress
(see below) and :math:`c` a constant offset. For example, for a cubic crystal the pattern
:math:`\mathbf e = (1, -1, 0, 0, 0, 0)` gives :math:`a = C_{11} - C_{12}`.

The unstrained cell is the same for every deformation, so its energy :math:`E_0` is calculated
once and used for the :math:`\delta = 0` point of every curve.

Applying the strain
-------------------

With the lattice vectors as rows of :math:`L`, the deformed lattice is :math:`L' = L\,F` with

.. math::

   F = I + \boldsymbol\varepsilon =
   \begin{pmatrix}
   1+\varepsilon_1 & \varepsilon_6/2 & \varepsilon_5/2\\
   \varepsilon_6/2 & 1+\varepsilon_2 & \varepsilon_4/2\\
   \varepsilon_5/2 & \varepsilon_4/2 & 1+\varepsilon_3
   \end{pmatrix},

the symmetric strain tensor built from the Voigt strain. The atoms keep their fractional
coordinates, so they move with the cell; with ``relaxation=True`` (the default; ``ISIF = 2`` in
VASP run folders) they are then relaxed inside the fixed strained cell, which gives the
relaxed-ion constants. ``relaxation=False`` gives the clamped-ion constants.

Volume-conserving corrections
-----------------------------

Six deformations get one extra diagonal term of second order in :math:`\delta`, chosen so that
the cell volume stays exactly the same, :math:`\det F = 1`. The term changes the energy only
at order :math:`\delta^3` and leaves the curvature :math:`a` unchanged; it removes the effect
of a small residual pressure :math:`P`, which would otherwise add
:math:`P\,\Delta V \propto \delta^2` to the curvature. The deformation names are those of
:doc:`guide/systems`.

**Opposite normal strains in x and y** – ``C11_C12_II`` (hexagonal, tetragonal, trigonal),
``C11_C22_C12_I`` (orthorhombic, monoclinic):

.. math::

   F = \begin{pmatrix} 1+\delta & 0 & 0 \\ 0 & 1-\delta & 0 \\ 0 & 0 & 1+\varepsilon_3 \end{pmatrix},
   \qquad \varepsilon_3 = \frac{\delta^2}{1-\delta^2}

**Opposite normal strains in x and z** – ``C11_C33_C13_I`` (hexagonal, tetragonal,
orthorhombic, monoclinic), ``C11_C33_C13_II`` (trigonal):

.. math::

   F = \begin{pmatrix} 1+\delta & 0 & 0 \\ 0 & 1+\varepsilon_2 & 0 \\ 0 & 0 & 1-\delta \end{pmatrix},
   \qquad \varepsilon_2 = \frac{\delta^2}{1-\delta^2}

**Opposite normal strains in y and z** – ``C22_C33_C23_I`` (orthorhombic, monoclinic):

.. math::

   F = \begin{pmatrix} 1+\varepsilon_1 & 0 & 0 \\ 0 & 1+\delta & 0 \\ 0 & 0 & 1-\delta \end{pmatrix},
   \qquad \varepsilon_1 = \frac{\delta^2}{1-\delta^2}

**Shear in the yz plane**, :math:`\varepsilon_4 = \delta` – ``C44`` (orthorhombic,
monoclinic):

.. math::

   F = \begin{pmatrix} 1+\varepsilon_1 & 0 & 0 \\ 0 & 1 & \delta/2 \\ 0 & \delta/2 & 1 \end{pmatrix},
   \qquad \varepsilon_1 = \frac{\delta^2}{4-\delta^2}

**Shear in the xz plane**, :math:`\varepsilon_5 = \delta` – ``C55`` (orthorhombic):

.. math::

   F = \begin{pmatrix} 1 & 0 & \delta/2 \\ 0 & 1+\varepsilon_2 & 0 \\ \delta/2 & 0 & 1 \end{pmatrix},
   \qquad \varepsilon_2 = \frac{\delta^2}{4-\delta^2}

**Shear in the xy plane**, :math:`\varepsilon_6 = \delta` – ``C66`` (tetragonal,
orthorhombic):

.. math::

   F = \begin{pmatrix} 1 & \delta/2 & 0 \\ \delta/2 & 1 & 0 \\ 0 & 0 & 1+\varepsilon_3 \end{pmatrix},
   \qquad \varepsilon_3 = \frac{\delta^2}{4-\delta^2}

In each case :math:`\det F = 1`: for the normal strains
:math:`(1+\delta)(1-\delta)(1+\varepsilon) = 1`, for the shears
:math:`(1+\varepsilon)(1-\delta^2/4) = 1`. The other deformations (all of the cubic table,
for example) have no such term.

Solving for the constants
-------------------------

The independent constants :math:`\mathbf p` of a table enter :math:`C` linearly,
:math:`C = \sum_k p_k B_k`, so every deformation gives one linear equation

.. math::

   a_j = \sum_k \left(\tfrac12\, \mathbf e_j^{\mathsf T} B_k\, \mathbf e_j\right) p_k ,

and the deformations of each table are chosen so that this system has full rank. It is solved
by least squares (:doc:`guide/systems` lists every :math:`a_j`).

Symmetry
--------

A rotation :math:`R` of the crystal transforms the Voigt stiffness matrix as
:math:`C' = M C M^{\mathsf T}`, with the 6×6 matrix

.. math::

   M_{IJ} = R_{ik} R_{jl} \;(J \le 3), \qquad
   M_{IJ} = R_{ik} R_{jl} + R_{il} R_{jk} \;(J > 3),

for the index pairs :math:`I \leftrightarrow (ij)`, :math:`J \leftrightarrow (kl)`. The
allowed stiffness matrices of a crystal are those with :math:`M C M^{\mathsf T} = C` for
every operation of its point group. This gives the orientation check (can a table represent
all of them?) and the orientation-aware automatic mode (a basis of them, and deformations
chosen to determine its coefficients). It gives 21, 13, 9, 7, 6, 7, 6, 5 and 3 independent
constants for the Laue classes -1, 2/m, mmm, 4/m, 4/mmm, -3, -3m, 6/mmm and m-3m.

Residual stress
---------------

For a reference cell under a residual stress :math:`\boldsymbol\sigma_0`, the energy gains a
linear term, :math:`(E - E_0)/V_0 = \boldsymbol\sigma_0\cdot\boldsymbol\varepsilon + \dots`, so
the linear coefficient of each fit is :math:`b_j = \boldsymbol\sigma_0 \cdot \mathbf e_j`.
Solving these equations for all deformations gives the components of
:math:`\boldsymbol\sigma_0` they determine; the residual pressure
:math:`P = -(\sigma_1 + \sigma_2 + \sigma_3)/3` is determined exactly when
:math:`(1,1,1,0,0,0)` is a combination of the :math:`\mathbf e_j`, which holds for every named
table. Because the cubic term of :math:`E(\delta)` leaks into the linear coefficient of a
quadratic fit, :math:`b_j` is taken from a cubic fit when there are at least four strains.

Mechanical stability
--------------------

A crystal is mechanically stable (Born criteria) when the strain energy is positive for
every small strain, i.e. :math:`C` is positive definite: all its eigenvalues are positive.
For cubic crystals this is :math:`C_{11} - C_{12} > 0`, :math:`C_{11} + 2C_{12} > 0`,
:math:`C_{44} > 0`; for hexagonal crystals :math:`C_{11} > |C_{12}|`,
:math:`2C_{13}^2 < C_{33}(C_{11}+C_{12})`, :math:`C_{44} > 0`, :math:`C_{66} > 0`.
elastic_ase decides with the eigenvalues, which is exact in any orientation.

Polycrystalline averages
------------------------

With the compliance :math:`S = C^{-1}`:

.. math::

   9K_V &= (C_{11}+C_{22}+C_{33}) + 2(C_{12}+C_{13}+C_{23})\\
   15G_V &= (C_{11}+C_{22}+C_{33}) - (C_{12}+C_{13}+C_{23}) + 3(C_{44}+C_{55}+C_{66})\\
   1/K_R &= (S_{11}+S_{22}+S_{33}) + 2(S_{12}+S_{13}+S_{23})\\
   15/G_R &= 4(S_{11}+S_{22}+S_{33}) - 4(S_{12}+S_{13}+S_{23}) + 3(S_{44}+S_{55}+S_{66})

The Hill values are the means :math:`K_H = (K_V+K_R)/2`, :math:`G_H = (G_V+G_R)/2`, and

.. math::

   E_H = \frac{9K_H G_H}{3K_H + G_H}, \qquad
   \nu_H = \frac{3K_H - 2G_H}{2(3K_H + G_H)}, \qquad
   A^U = 5\frac{G_V}{G_R} + \frac{K_V}{K_R} - 6 .

Sound velocities and Debye temperature
--------------------------------------

For the isotropic polycrystal with density :math:`\rho`:

.. math::

   v_l = \sqrt{\frac{K_H + \frac43 G_H}{\rho}}, \qquad
   v_t = \sqrt{\frac{G_H}{\rho}}, \qquad
   v_m = \left[\frac13\left(\frac{2}{v_t^3} + \frac{1}{v_l^3}\right)\right]^{-1/3},

and the Debye temperature (Anderson, J. Phys. Chem. Solids 24, 909 (1963))

.. math::

   \theta_D = \frac{h}{k_B}\left(\frac{3n}{4\pi V}\right)^{1/3} v_m ,

with :math:`n/V` the number of atoms per volume.

Hardness estimates
------------------

With :math:`k = G_H/K_H`: Chen et al., Intermetallics 19, 1275 (2011),
:math:`H_V = 2(k^2 G_H)^{0.585} - 3`; Tian et al., Int. J. Refract. Met. Hard Mater. 33, 93
(2012), :math:`H_V = 0.92\,k^{1.137} G_H^{0.708}` (GPa). Both are empirical fits to hard
covalent and ionic solids.

Limitations
-----------

* The reference structure must be relaxed with the same calculator and settings (the
  residual pressure is reported).
* Anharmonic terms bias the quadratic fit at larger strains (a few tenths of a percent at
  ±3 % for typical materials).
* The symmetry analysis ignores magnetic order.
* The constants are static values of the calculator's energy surface; zero-point and thermal
  effects are not included.
