"""
Command-line interface (CLI) for ElasticASE.
Allows calculating elastic constants directly from terminal for any structure file.
"""

import os
import sys
import argparse
import json
import numpy as np
import pandas as pd
from ase.io import read
from elastic_ase.calculator import ElasticASE

def parse_args(args=None):
    parser = argparse.ArgumentParser(
        prog="elastic-ase",
        description="Compute second-order elastic constants (C_ij) from structure files with ASE.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "structure",
        type=str,
        help="Path to atomic structure file (CIF, POSCAR, XYZ, etc. readable by ASE)."
    )
    parser.add_argument(
        "-c", "--calculator",
        type=str,
        default="emt",
        choices=["emt", "grace"],
        help="ASE calculator to use ('emt' for built-in test potential, 'grace' for GRACE foundation models)."
    )
    parser.add_argument(
        "-m", "--model",
        type=str,
        default="GRACE-2L-OAM",
        help="GRACE model name if calculator is 'grace'."
    )
    parser.add_argument(
        "-s", "--symmetry",
        type=str,
        default="auto",
        choices=["auto", "cubic", "hexagonal", "orthogonal", "trigonal_1", "trigonal_2", "monoclinic", "triclinic"],
        help="Crystal symmetry family."
    )
    parser.add_argument(
        "--strains",
        type=str,
        default="-0.01,-0.005,0.0,0.005,0.01",
        help="Comma-separated strain values."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Optional output file path to save results (JSON or CSV format based on extension)."
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress detailed terminal output."
    )
    return parser.parse_args(args)

def main(args=None):
    parsed = parse_args(args)

    if not os.path.exists(parsed.structure):
        print(f"Error: Structure file '{parsed.structure}' not found.", file=sys.stderr)
        return 1

    # Load structure
    atoms = read(parsed.structure)

    # Attach calculator
    if parsed.calculator == "emt":
        from ase.calculators.emt import EMT
        atoms.calc = EMT()
    elif parsed.calculator == "grace":
        os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
        from tensorpotential.calculator import grace_fm
        atoms.calc = grace_fm(parsed.model)

    # Parse strains
    strain_list = [float(x.strip()) for x in parsed.strains.split(",") if x.strip()]

    # Run ElasticASE
    el = ElasticASE(
        atoms,
        symmetry=parsed.symmetry,
        strains=strain_list,
        verbose=not parsed.quiet
    )
    results = el.calculate()

    # Save output if requested
    if parsed.output:
        out_path = parsed.output
        if out_path.endswith(".json"):
            serializable = {
                "symmetry": str(el.symmetry),
                "stiffness_matrix": np.array(el.stiffness).tolist(),
                "is_stable": bool(el.is_stable),
                "stability_criteria": str(el.stability_criteria),
                "polycrystalline_averages": {k: float(v) for k, v in el.averages.items()},
                "r2_scores": {k: float(v) for k, v in el.r2_scores.items()}
            }
            with open(out_path, "w") as f:
                json.dump(serializable, f, indent=2)
            if not parsed.quiet:
                print(f"\n[ElasticASE] Saved results to {out_path}")
        elif out_path.endswith(".csv"):
            data = {
                "Symmetry": [el.symmetry],
                "Stable": [el.is_stable],
                "K_VRH_GPa": [el.averages["K_VRH"]],
                "G_VRH_GPa": [el.averages["G_VRH"]],
                "E_VRH_GPa": [el.averages["E_VRH"]],
                "nu": [el.averages["nu_VRH"]],
                "Pugh_ratio": [el.averages["Pugh_ratio"]]
            }
            for i in range(6):
                for j in range(6):
                    data[f"C{i+1}{j+1}_GPa"] = [round(float(el.stiffness[i, j]), 2)]
            pd.DataFrame(data).to_csv(out_path, index=False)
            if not parsed.quiet:
                print(f"\n[ElasticASE] Saved results to {out_path}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
