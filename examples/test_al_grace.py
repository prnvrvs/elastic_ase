import os
import sys

# Suppress TF logs
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

from ase.build import bulk
from tensorpotential.calculator import grace_fm
from elastic_ase import ElasticASE

def main():
    print("Testing elastic_ase on FCC Aluminium with GRACE-2L-OAM...")
    calc = grace_fm("GRACE-2L-OAM")
    al = bulk("Al", "fcc", a=4.0428, cubic=True)
    al.calc = calc

    el = ElasticASE(al, symmetry="cubic")
    res = el.calculate()

    return res

if __name__ == "__main__":
    main()
