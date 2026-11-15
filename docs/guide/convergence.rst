Convergence tests
=================

Elastic constants come from very small energy differences. For fcc aluminium at ±1 % strain
the energy changes by about 0.2–2 meV per atom, depending on the deformation, while the total energy is usually converged
to about 1 meV per atom. Settings that are "converged" for the total energy can therefore still
give wrong elastic constants. **Converge the elastic constants themselves.**

What to converge, in order of importance:

* **k-point mesh** – the most important for metals, whose shear constants can need very
  dense meshes;
* **Brillouin-zone integration** – the smearing method and width (VASP ``ISMEAR`` /
  ``SIGMA``). The tetrahedron method with Blöchl corrections (``ISMEAR = -5``) gives
  inaccurate forces: do not use it when the atoms are relaxed;
* **plane-wave cutoff** – usually much less critical than the k-mesh;
* **SCF tolerance** – tight enough that the energy noise is far below the strain energy
  (e.g. VASP ``EDIFF = 1E-8``);
* **strain range** – the default ±1 % suits most calculations; see
  `Choosing the strain range`_ below.

Relax the structure with the final settings, and check the residual pressure that elastic_ase
reports and the energy–strain curves (``plot_fits``).

Doing it with elastic_ase
-------------------------

Run the same structure for a series of k-meshes, each in its own work folder, and compare the
constants (VASP run folders; the same works with an ASE calculator and ``calculate()``):

.. code-block:: python

   from ase.io import read
   from elastic_ase import ElasticConstants

   atoms = read("CONTCAR")                    # relaxed structure
   results = {}
   for k in (12, 16, 20, 24, 28, 32):
       kpoints = f"KPOINTS_{k}"
       with open(kpoints, "w") as f:
           f.write(f"{k}x{k}x{k}\n0\nGamma\n {k} {k} {k}\n 0 0 0\n")
       ec = ElasticConstants(atoms, code="VASP",
                             command="mpirun -np 8 vasp_std", verbose=False)
       work = f"kconv/k{k}"
       files = {"INCAR": "INCAR", "KPOINTS": kpoints, "POTCAR": "POTCAR"}
       ec.write_inputs(work, files=files)
       ec.run_jobs(work)                      # or submit the folders yourself
       results[k] = ec.read_results(work)

   for k, r in results.items():
       c = r.independent_constants
       print(k, {name: round(v, 1) for name, v in c.items()},
             "K = %.1f GPa" % r.polycrystalline_averages["K_VRH"],
             "E0 = %.5f eV" % r.energies["ground_state"])

Plot the constants against the mesh and choose the mesh beyond which they change by less than
the accuracy you need (for example 1–2 GPa):

.. code-block:: python

   import matplotlib.pyplot as plt

   ks = sorted(results)
   for name in results[ks[0]].independent_constants:
       values = [results[k].independent_constants[name] for k in ks]
       plt.plot(ks, values, "o-", label=name)
   plt.xlabel("k-mesh (n x n x n)")
   plt.ylabel("GPa")
   plt.legend()
   plt.savefig("k_convergence.png", dpi=150)

Repeat for the cutoff (``ENCUT``), the smearing and the strain range (``strains=[...]``) once
the k-mesh is fixed. Each run is resumable: calling the loop again only repeats missing runs.

.. _choosing-the-strain-range:

Choosing the strain range
-------------------------

The strain range trades two errors against each other:

* **Noise.** The strain energy grows with :math:`\delta^2`: at ±2 % the energy differences
  are 4 times larger than at ±1 %, so numerical noise of the calculator matters 4 times less.
* **Anharmonic bias.** A parabola is fitted to the energies. The higher-order terms of the
  energy also grow with the strain range and shift each constant steadily in one direction as
  the range grows.

The bias, for GRACE-2L-OAM with the atoms relaxed in every strained cell (5 strains from
:math:`-\delta` to :math:`+\delta`; constants in GPa):

.. list-table::
   :header-rows: 1
   :widths: 22 13 13 13 13 13 13

   * - Constant
     - ±0.5 %
     - ±1 %
     - ±1.5 %
     - ±2 %
     - ±3 %
     - Change, ±1 % → ±3 %
   * - Al :math:`C_{11}`
     - 112.80
     - 112.84
     - 112.90
     - 112.99
     - 113.22
     - +0.3 %
   * - Al :math:`C_{12}`
     - 60.48
     - 60.53
     - 60.61
     - 60.71
     - 61.01
     - +0.8 %
   * - Al :math:`C_{44}`
     - 33.86
     - 33.86
     - 33.87
     - 33.89
     - 33.93
     - +0.2 %
   * - Si :math:`C_{11}`
     - 143.22
     - 142.76
     - 142.31
     - 141.94
     - 141.36
     - −1.0 %
   * - Si :math:`C_{12}`
     - 68.09
     - 67.70
     - 67.37
     - 67.17
     - 66.99
     - −1.0 %
   * - Si :math:`C_{44}`
     - 76.35
     - 76.24
     - 76.04
     - 75.80
     - 75.25
     - −1.3 %
   * - Mg :math:`C_{11}`
     - 61.86
     - 61.87
     - 61.88
     - 61.89
     - 61.94
     - +0.1 %
   * - Mg :math:`C_{12}`
     - 27.24
     - 27.25
     - 27.27
     - 27.30
     - 27.38
     - +0.5 %
   * - Mg :math:`C_{13}`
     - 21.84
     - 21.83
     - 21.82
     - 21.84
     - 21.91
     - +0.4 %
   * - Mg :math:`C_{33}`
     - 65.81
     - 65.83
     - 65.88
     - 66.03
     - 66.48
     - +1.0 %
   * - Mg :math:`C_{44}`
     - 18.11
     - 18.10
     - 18.07
     - 18.03
     - 17.94
     - −0.9 %
   * - TiO₂ :math:`C_{11}`
     - 261.91
     - 262.04
     - 262.23
     - 262.47
     - 263.03
     - +0.4 %
   * - TiO₂ :math:`C_{12}`
     - 170.73
     - 170.57
     - 170.32
     - 170.01
     - 169.28
     - −0.8 %
   * - TiO₂ :math:`C_{13}`
     - 143.81
     - 144.07
     - 144.52
     - 145.15
     - 146.88
     - +2.0 %
   * - TiO₂ :math:`C_{33}`
     - 459.93
     - 460.32
     - 461.21
     - 462.72
     - 467.25
     - +1.5 %
   * - TiO₂ :math:`C_{44}`
     - 103.48
     - 103.49
     - 103.51
     - 103.54
     - 103.64
     - +0.1 %
   * - TiO₂ :math:`C_{66}`
     - 210.14
     - 210.17
     - 210.22
     - 210.29
     - 210.49
     - +0.2 %

What to use:

* **±1 % (default)** for most calculations. With a smooth calculator (machine-learning or
  classical potential) ±0.5 % is equally good.
* **±2 %** when the energies are noisy – a coarse k-mesh for a metal, a loose SCF tolerance.
  The bias stays below about 1 % of the constants.
* **Unsure?** Run two ranges, e.g. ``strains=[-0.01, -0.005, 0, 0.005, 0.01]`` and
  ``strains=[-0.02, -0.01, 0, 0.01, 0.02]``. If the constants agree to the accuracy you need,
  either is fine; if not, the smaller range is closer when its curves are smooth.

