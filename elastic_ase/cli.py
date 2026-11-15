"""Command line tool ``elastic-ase``: elastic constants of a structure file."""

import os
import sys
import argparse
import csv
import importlib
import importlib.util
import json
import numpy as np
from ase.io import read
from elastic_ase.calculator import ElasticASE, _plain
from elastic_ase.symmetry import check_symmetry


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        prog="elastic-ase",
        description="Compute second-order elastic constants (C_ij) and visualize anisotropic properties with ASE.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "structure",
        type=str,
        help="Path to atomic structure file (CIF, POSCAR, XYZ, etc. readable by ASE)."
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Perform a crystal symmetry and orientation check without calculating elastic constants."
    )
    parser.add_argument(
        "-c", "--calculator",
        type=str,
        default="emt",
        help="ASE calculator: 'emt' (built-in test potential), 'grace' (GRACE foundation model, "
             "see --model), or any calculator as 'module:function', e.g. "
             "'mace.calculators:mace_mp' (called without arguments)."
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
        choices=[
            "auto", "cubic", "hexagonal", "orthogonal", "orthorhombic",
            "tetragonal_1", "tetragonal_2", "trigonal_1", "trigonal_2",
            "monoclinic", "triclinic"
        ],
        help="Deformation table: 'auto' (default) or a named table."
    )
    parser.add_argument(
        "--frame",
        type=str,
        default="cartesian",
        choices=["cartesian", "standard"],
        help="Coordinate frame for output stiffness tensor ('cartesian' or 'standard')."
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
        help="Save the results: .json (everything), .csv, any other extension: plain 6x6 matrix."
    )
    parser.add_argument(
        "--no-reorient",
        dest="reorient",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Keep the structure in its given orientation (default: rotate it rigidly into the "
             "standard setting; cell and atoms are never changed)."
    )
    parser.add_argument(
        "--no-relaxation",
        dest="relaxation",
        action="store_false",
        default=argparse.SUPPRESS,
        help="Keep the atoms at their scaled positions in the strained cells (clamped-ion "
             "constants). Default: relax them (relaxed-ion constants)."
    )
    parser.add_argument(
        "--fmax",
        type=float,
        default=1e-3,
        help="Force criterion (eV/Å) of the atom relaxation in the strained cells and the "
             "unstrained cell."
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="No report."
    )

    # Visualization flags powered by elastic-vis
    vis_group = parser.add_argument_group("Visualization options (elastic-vis)")
    vis_group.add_argument(
        "--plot", action="store_true",
        help="Generate directional elastic property plots via elastic-vis."
    )
    vis_group.add_argument("--Young", action="store_true", help="Young's modulus plots (default: Young and shear if none is chosen).")
    vis_group.add_argument("--Shear", action="store_true", help="Shear modulus plots.")
    vis_group.add_argument("--Poisson", action="store_true", help="Poisson's ratio plots.")
    vis_group.add_argument("--LC", action="store_true", help="Linear compressibility plots.")
    vis_group.add_argument("--2D", dest="plot_2d", action="store_true", help="2D polar plots (default: 2D and 3D if neither is chosen).")
    vis_group.add_argument("--3D", dest="plot_3d", action="store_true", help="3D surface plots.")
    vis_group.add_argument(
        "--format", type=str, default="png", choices=["png", "pdf", "svg"],
        help="Image format for saved figures."
    )
    vis_group.add_argument(
        "--output-prefix", type=str, default=None,
        help="Prefix for generated visualization files (defaults to structure file stem)."
    )
    vis_group.add_argument(
        "--style", type=str, default=None,
        help="Path to custom JSON styling file for elastic-vis."
    )
    vis_group.add_argument(
        "--show", action="store_true",
        help="Show plots interactively on screen."
    )

    args = list(sys.argv[1:] if args is None else args)
    for i, arg in enumerate(args[:-1]):
        if arg == "--strains" and args[i + 1].startswith("-"):
            args[i:i + 2] = [f"--strains={args[i + 1]}"]
            break
    return parser.parse_args(args)


def make_calculator(name, model=None):
    """ASE calculator from a CLI name: "emt", "grace" or "module:function"."""
    if name == "emt":
        from ase.calculators.emt import EMT
        return EMT()
    if name == "grace":
        os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
        from tensorpotential.calculator import grace_fm
        return grace_fm(model)
    if ":" in name:
        module_name, function_name = name.split(":", 1)
        return getattr(importlib.import_module(module_name), function_name)()
    raise ValueError("use 'emt', 'grace' or 'module:function'")


def main(args=None):
    parsed = parse_args(args)

    if not os.path.exists(parsed.structure):
        print(f"Error: Structure file '{parsed.structure}' not found.", file=sys.stderr)
        return 1

    # Fast symmetry check mode
    if parsed.check:
        check_symmetry(parsed.structure, verbose=True)
        return 0

    # Load structure
    atoms = read(parsed.structure)

    # Attach calculator
    try:
        atoms.calc = make_calculator(parsed.calculator, parsed.model)
    except (ImportError, AttributeError, ValueError, TypeError) as exc:
        hint = " (a 'module:function' calculator must work without arguments)" if isinstance(exc, TypeError) else ""
        print(f"Error: cannot create calculator '{parsed.calculator}': {exc}{hint}", file=sys.stderr)
        return 1

    if parsed.plot and importlib.util.find_spec("elastic_vis") is None:
        print("Error: --plot needs elastic-vis: "
              'pip install "git+https://github.com/prabhat7iitm/elastic-vis.git"', file=sys.stderr)
        return 1

    # Parse strains and set up
    try:
        strain_list = [float(x.strip()) for x in parsed.strains.split(",") if x.strip()]
        el = ElasticASE(
            atoms,
            symmetry=parsed.symmetry,
            strains=strain_list,
            reorient=getattr(parsed, "reorient", True),
            relaxation=getattr(parsed, "relaxation", True),
            fmax=parsed.fmax,
            frame=parsed.frame,
            verbose=not parsed.quiet
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    results = el.calculate()

    # Save output if requested
    if parsed.output:
        out_path = parsed.output
        active_matrix = results["stiffness_matrix"]
        if out_path.endswith(".json"):
            serializable = {
                "symmetry": str(el.symmetry),
                "frame": str(el.frame),
                "stiffness_matrix": np.array(active_matrix).tolist(),
                "stiffness_cartesian": np.array(el.stiffness).tolist() if el.stiffness is not None else None,
                "stiffness_standard": np.array(el.stiffness_standard).tolist() if el.stiffness_standard is not None else None,
                "standard_constants": {k: float(v) for k, v in el.standard_constants.items()} if el.standard_constants else {},
                "compliance_matrix": np.array(el.compliance).tolist(),
                "is_stable": bool(el.is_stable),
                "stability_criteria": str(el.stability_criteria),
                "polycrystalline_averages": {k: float(v) for k, v in el.averages.items()},
                "directional_properties": {k: float(v) if isinstance(v, (int, float, np.floating)) else v for k, v in el.directional.items()} if el.directional else None,
                "r2_scores": {k: float(v) for k, v in el.r2_scores.items()},
                "independent_constants": _plain(results["independent_constants"]),
                "derived_properties": _plain(results["derived_properties"]),
                "residual_stress": _plain(results["residual_stress"]),
                "energies_eV": _plain(results["energies"]),
            }
            with open(out_path, "w") as f:
                json.dump(serializable, f, indent=2)
            if not parsed.quiet:
                print(f"\n[ElasticASE] Saved results to {out_path}")
        elif out_path.endswith(".csv"):
            data = {
                "Symmetry": [el.symmetry],
                "Frame": [el.frame],
                "Stable": [el.is_stable],
                "K_VRH_GPa": [el.averages["K_VRH"]],
                "G_VRH_GPa": [el.averages["G_VRH"]],
                "E_VRH_GPa": [el.averages["E_VRH"]],
                "nu": [el.averages["nu_VRH"]],
                "Pugh_ratio": [el.averages["Pugh_ratio"]],
                "A_universal": [el.averages["universal_anisotropy"]]
            }
            if el.directional:
                data["E_min_GPa"] = [el.directional["E_min"]]
                data["E_max_GPa"] = [el.directional["E_max"]]
                data["E_anisotropy"] = [el.directional["E_anisotropy"]]
                data["G_min_GPa"] = [el.directional["G_min"]]
                data["G_max_GPa"] = [el.directional["G_max"]]
                data["G_anisotropy"] = [el.directional["G_anisotropy"]]
                data["nu_min"] = [el.directional["nu_min"]]
                data["nu_max"] = [el.directional["nu_max"]]
                data["Auxetic"] = [el.directional["is_auxetic"]]
            for i in range(6):
                for j in range(6):
                    data[f"C{i+1}{j+1}_GPa"] = [round(float(active_matrix[i, j]), 2)]
            with open(out_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(list(data))
                writer.writerow([column[0] for column in data.values()])
            if not parsed.quiet:
                print(f"\n[ElasticASE] Saved results to {out_path}")
        else:
            # Default to plain 6x6 text matrix (readable by elastic-vis & numpy.loadtxt)
            el.export_matrix(out_path)

    # Visualizations via elastic-vis
    if parsed.plot:
        props = []
        if parsed.Young: props.append("Young")
        if parsed.Shear: props.append("Shear")
        if parsed.Poisson: props.append("Poisson")
        if parsed.LC: props.append("LC")
        if not props: props = ["Young", "Shear"]

        dims = []
        if parsed.plot_2d: dims.append("2D")
        if parsed.plot_3d: dims.append("3D")
        if not dims: dims = ["2D", "3D"]

        prefix = parsed.output_prefix or os.path.splitext(os.path.basename(parsed.structure))[0]
        el.visualize(
            properties=props,
            dimensions=dims,
            prefix=prefix,
            style_config=parsed.style,
            show=parsed.show,
            output_format=parsed.format
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
