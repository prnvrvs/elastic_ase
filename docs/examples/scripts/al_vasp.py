"""Elastic constants of fcc aluminium with VASP, in run folders.

Needs: pip install elastic-ase; VASP (vasp_std); and the PBE POTCAR of Al in this folder
Run:   python al_vasp.py
The settings are a starting point; converge ENCUT and the k-point mesh on the elastic constants.
"""
from ase.build import bulk

from elastic_ase import ElasticConstants

# Step 1: settings and structure
atoms = bulk("Al", "fcc", a=4.04)   # use your relaxed structure, e.g. read("CONTCAR")

with open("INCAR", "w") as f:
    f.write("""ENCUT  = 520
EDIFF  = 1E-8
ISMEAR = 1
SIGMA  = 0.1
IBRION = 2
ISIF   = 2
NSW    = 50
EDIFFG = -1E-3
LWAVE  = .FALSE.
LCHARG = .FALSE.
""")
with open("KPOINTS", "w") as f:
    f.write("Gamma-centred mesh\n0\nGamma\n 24 24 24\n 0 0 0\n")

# Step 2: write, run and read
ec = ElasticConstants(atoms, code="VASP", command="mpirun -np 16 vasp_std")
ec.write_inputs("runs", files={"INCAR": "INCAR", "KPOINTS": "KPOINTS", "POTCAR": "POTCAR"})
ec.run_jobs("runs")                 # or submit the folders on your cluster
ec.vasp_status("runs")              # done / failed / pending
results = ec.read_results("runs")
ec.save("al_vasp.json")
