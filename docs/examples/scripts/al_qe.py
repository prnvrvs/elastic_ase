"""Elastic constants of fcc aluminium with Quantum ESPRESSO (pw.x), in run folders.

Needs: pip install elastic-ase; pw.x (e.g. conda install -c conda-forge qe); and the
       pseudopotential Al.pbe-n-kjpaw_psl.1.0.0.UPF in ./pseudo (from pslibrary or the
       Quantum ESPRESSO pseudopotential page)
Run:   python al_qe.py
The settings are a quick start (12^3 k-points), not converged for Al.
"""
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile

from elastic_ase import ElasticConstants

atoms = bulk("Al", "fcc", a=4.04)
atoms.calc = Espresso(
    profile=EspressoProfile(command="pw.x", pseudo_dir="pseudo"),
    pseudopotentials={"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
    input_data={"control": {"tprnfor": True, "tstress": True},
                "system": {"ecutwfc": 40, "ecutrho": 320, "occupations": "smearing",
                           "smearing": "mv", "degauss": 0.02},
                "electrons": {"conv_thr": 1e-10}},
    kpts=(12, 12, 12))

ec = ElasticConstants(atoms, code="qe",
                      command="mpirun -np 6 pw.x -in espresso.pwi > espresso.pwo")
ec.write_inputs("runs")
ec.run_jobs("runs")
print(ec.qe_status("runs"))
results = ec.read_results("runs")
ec.save("al_qe.json")
