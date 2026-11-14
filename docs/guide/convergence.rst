Convergence tests
=================

Elastic constants come from very small energy differences. For fcc aluminium at ±1 % strain
the energy changes by less than 1 meV per atom, while the total energy is usually converged
to about 1 meV per atom. Settings that are "converged" for the total energy can therefore still
give wrong elastic constants. **Converge the elastic constants themselves.**

What to converge, in order of importance:

* **k-point mesh** – the most important for metals, whose shear constants can need very
  dense meshes;
* **Brillouin-zone integration** – the smearing method and width (VASP ``ISMEAR`` /
  ``SIGMA``); the tetrahedron method usually converges faster than smearing for elastic
  constants of metals;
* **plane-wave cutoff** – usually much less critical than the k-mesh;
* **SCF tolerance** – tight enough that the energy noise is far below the strain energy
  (e.g. VASP ``EDIFF = 1E-8``);
* **strain range** – the default ±1 % keeps anharmonic terms small; larger strains raise the
  signal (the strain energy grows with :math:`\delta^2`) but bring in anharmonic terms. Compare
  both.

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
       ec = ElasticConstants(atoms, code="vasp",
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
