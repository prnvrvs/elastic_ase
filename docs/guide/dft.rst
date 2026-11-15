DFT codes: VASP, QE, …
======================

A DFT calculation takes minutes to hours, can fail, and usually runs on a cluster. elastic_ase
therefore works with **run folders**: one folder per calculation, written by elastic_ase, run
by the DFT code (by ``run_jobs`` or by your own job scripts), and read back by elastic_ase.

Tell ``ElasticConstants`` once which code you use; every later step then knows what to write,
which output file to read and check, and which command to run.

.. code-block:: python

   from ase.io import read
   from elastic_ase import ElasticConstants

   ec = ElasticConstants(read("POSCAR"), code="VASP",
                         command="mpirun -np 16 vasp_std")
   files = {"INCAR": "INCAR", "KPOINTS": "KPOINTS", "POTCAR": "POTCAR"}
   ec.write_inputs("runs", files=files)
   ec.run_jobs("runs")
   results = ec.read_results("runs")

.. list-table:: What ``code=`` sets
   :header-rows: 1
   :widths: 12 34 32 22

   * - ``code``
     - ``write_inputs`` writes
     - ``status`` / ``read_results`` read and check
     - default ``command``
   * - ``"VASP"``
     - ``POSCAR`` + your ``files`` (or the full input from an ASE ``Vasp`` calculator)
     - ``OUTCAR``: readable, expected cell, electronic loop converged
     - ``vasp_std``
   * - ``"QE"``
     - ``espresso.pwi`` from an ASE ``Espresso`` calculator
     - ``espresso.pwo``: readable, expected cell, SCF converged
     - ``pw.x -in espresso.pwi > espresso.pwo``
   * - none
     - the attached calculator's own input if it is a file-based ASE calculator, otherwise
       the structure (``format=``); plus your ``files``
     - the file given as ``output=`` (anything ``ase.io.read`` reads)
     - the one given as ``command=``

``output="vasprun.xml"`` reads a file other than the code's default; the VASP convergence
check works only on OUTCAR.

The run folders
---------------

.. code-block:: text

   runs/
     ground_state/                 unstrained cell, run once
     C11_C12_I/strain_-0.010000/   one folder per deformation and strain
     C11_C12_I/strain_-0.005000/
     ...
     elastic_ase_runs.json         symmetry, strains and list of runs

Run names such as ``C44/strain_+0.005000`` are what ``status``, ``missing`` and the ``runs=``
arguments use.

VASP
----

1. Relax the structure with your converged settings (``ISIF = 3``), and start from its
   ``CONTCAR``.
2. Prepare one INCAR for all runs. In a multi-atom cell the atoms move when the cell is
   strained, so relax them inside the fixed strained cell (``ISIF = 2``) to get the
   relaxed-ion constants:

   .. code-block:: text

      ENCUT  = 520        ! your converged values
      EDIFF  = 1E-8       ! tight: the energy differences are small
      ISMEAR = 1
      SIGMA  = 0.1
      IBRION = 2          ! relax the atoms
      ISIF   = 2          ! keep the strained cell fixed
      NSW    = 50         ! maximum relaxation steps
      EDIFFG = -1E-3      ! force criterion (eV/Å)
      LWAVE  = .FALSE.
      LCHARG = .FALSE.

   Use the same KPOINTS and POTCAR for every run. In run
   folders the ``relaxation`` option of ``ElasticConstants`` has no effect: the INCAR
   decides.
3. Write, run, read:

   .. code-block:: python

      ec = ElasticConstants(read("CONTCAR"), code="VASP",
                            command="mpirun -np 16 vasp_std")
      files = {"INCAR": "INCAR", "KPOINTS": "KPOINTS", "POTCAR": "POTCAR"}
      ec.write_inputs("runs", files=files)
      ec.run_jobs("runs")        # or submit the folders on your cluster
      ec.status("runs")
      results = ec.read_results("runs")
      ec.plot_fits("energy_strain.png")

If VASP stops with "I REFUSE TO CONTINUE WITH THIS SICK JOB" because its symmetry analysis
finds the direct and reciprocal lattices of a strained cell inconsistent – this can happen for
cells kept in a non-standard orientation (``reorient=False``) – add ``ISYM = 0`` (or change
``SYMPREC``) to the INCAR and rerun the failed runs. The default rotation into the standard
setting usually avoids it.

The POTCAR's element order is checked against the POSCAR, and a structure whose atoms are not
grouped by element (e.g. Cu Au Cu) is refused, since no POTCAR can match its POSCAR – sort it
first with ``atoms = atoms[atoms.numbers.argsort(kind="stable")]``.

With an ASE ``Vasp`` calculator attached and no ``files``, ``write_inputs`` writes the complete
input (INCAR, KPOINTS, POSCAR, POTCAR) from the calculator instead (ASE needs
``VASP_PP_PATH`` pointing to your POTCAR library).

Shortcuts with the code fixed: ``vasp_input(workdir, incar=, kpoints=, potcar=)``,
``vasp_status(workdir)`` and ``vasp_output(workdir)``.

QE
--

You need a pw.x executable (for example ``conda install -c conda-forge qe``) and the
pseudopotentials, e.g. from the `SSSP library
<https://www.materialscloud.org/discover/sssp>`_ or pslibrary.

The input is written from an ASE ``Espresso`` calculator, whose pseudopotentials, k-points and
``input_data`` are used for every run:

.. code-block:: python

   from ase.calculators.espresso import Espresso, EspressoProfile

   atoms.calc = Espresso(
       profile=EspressoProfile(command="pw.x", pseudo_dir="pseudo"),
       pseudopotentials={"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
       input_data={"control": {"tprnfor": True, "tstress": True},
                   "system": {"ecutwfc": 60, "occupations": "smearing",
                              "smearing": "mv", "degauss": 0.02},
                   "electrons": {"conv_thr": 1e-10}},
       kpts=(16, 16, 16))

   cmd = "mpirun -np 16 pw.x -in espresso.pwi > espresso.pwo"
   ec = ElasticConstants(atoms, code="QE", command=cmd)
   ec.write_inputs("runs")
   ec.run_jobs("runs")
   results = ec.read_results("runs")

This is a single point in every run. To relax the atoms inside every fixed strained cell, add
``"calculation": "relax"`` (and, e.g., ``"forc_conv_thr": 1e-4``) to ``"control"``; elastic_ase reads the final energy of each run. As for VASP, the
``relaxation`` option of ``ElasticConstants`` has no effect in run folders.

A pw.x output containing ``convergence NOT achieved`` is reported as failed. Shortcuts:
``qe_input``, ``qe_status``, ``qe_output``.

Other codes
-----------

For a code without a ``code=`` setting, give the structure format, the extra input files, the
output file and the command yourself. This works for any code whose structure file ASE can
write and whose energy ASE can read from the output. Outputs other than VASP OUTCARs and pw.x
outputs get no convergence check:

.. code-block:: python

   ec = ElasticConstants(atoms, output="<output file>",
                         command="<command>")
   ec.write_inputs("runs", format="<ase format>",
                   filename="<structure file>",
                   files={"<input file>": "<path>"})
   ec.run_jobs("runs")
   results = ec.read_results("runs")

If the attached calculator is a file-based ASE calculator, its own input is written instead of
the structure; pass ``calculator_inputs=False`` to write only the structure.

Running, checking, repeating
----------------------------

``run_jobs(workdir)`` runs ``command`` with the shell inside every folder that has no usable
output yet, one after the other; the command's screen output goes to
``<folder>/elastic_ase_run.log`` and the exit codes are kept in ``ec.exit_codes``. It returns
``status(workdir)``.

.. code-block:: python

   ec.status("runs")       # done / failed (with the reason) / pending
   ec.is_complete("runs")  # every run has a usable output
   ec.missing("runs")      # failed + pending
   ec.run_jobs("runs")     # again: only the missing runs

A run is *failed* when its output exists but is unusable – unreadable, a different cell than
the folder should hold (mixed-up folders), or not converged – and *pending* when there is no
output yet.

On a cluster
------------

``run_jobs`` runs the folders one after the other. For many expensive runs, submit every
folder as its own job instead, for example with a job array over the folder paths that
``write_inputs`` returns:

.. code-block:: python

   folders = ec.write_inputs("runs", files={...})

and check and finish from Python afterwards – with the same structure, symmetry and strains.
``missing`` gives the runs to resubmit, and ``write_inputs(..., runs=...)`` rewrites only
those folders:

.. code-block:: python

   ec = ElasticConstants(read("CONTCAR"), code="VASP")
   ec.status("runs")
   ec.write_inputs("runs", files={...}, runs=ec.missing("runs"))
   # ... resubmit them ...
   results = ec.read_results("runs")

Use the same structure, ``symmetry`` and ``strains`` for writing and reading: the run names and
expected cells come from them, and a mismatch is reported.

Running DFT through ASE directly
--------------------------------

Instead of run folders, an ASE DFT calculator can also be driven by
``calculate(workdir=...)`` (:doc:`calculators`): every calculation runs in its own folder
inside ``workdir`` and every energy is saved as soon as it is computed, so a crashed job
continues where it stopped.

Convergence
-----------

Converge the elastic constants, not only the total energy, and use the same settings for every
run: :doc:`convergence`.
