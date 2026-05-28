"""
Generate conformers and Rosetta .params files from an SDF of ligands.

Usage:
    python make_conformers_and_params.py --input <file.sdf> --license-key <key>
                                         [--conformator-path <path>]
                                         [--num-conformers <n>] [--keep-3d]
                                         [--output-dir <dir>]
"""

import argparse
import os
import shutil
import sys
from pathlib import Path


CONFORMATOR_DEFAULT = "/conformator_for_container/conformator_1.2.1/conformator"
MOLFILE_TO_PARAMS = "/conformator_for_container/molfile_to_params.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate conformers and Rosetta params from an SDF file."
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to the input SDF file (up to ~5k ligands)."
    )
    parser.add_argument(
        "--license-key", "-k", required=True,
        help="Conformator license key."
    )
    parser.add_argument(
        "--conformator-path",
        default=CONFORMATOR_DEFAULT,
        help=f"Path to the conformator executable (default: {CONFORMATOR_DEFAULT})."
    )
    parser.add_argument(
        "--molfile-to-params",
        default=MOLFILE_TO_PARAMS,
        help=f"Path to molfile_to_params.py (default: {MOLFILE_TO_PARAMS})."
    )
    parser.add_argument(
        "--num-conformers", "-n", type=int, default=15,
        help="Max conformers per ligand (default: 15)."
    )
    parser.add_argument(
        "--keep-3d", action="store_true", default=True,
        help="Keep existing 3D coordinates (default)."
    )
    parser.add_argument(
        "--no-keep-3d", dest="keep_3d", action="store_false",
        help="Do not keep existing 3D coordinates."
    )
    parser.add_argument(
        "--hydrogens", action="store_true", default=True,
        help="Add hydrogens (default)."
    )
    parser.add_argument(
        "--no-hydrogens", dest="hydrogens", action="store_false",
        help="Do not add hydrogens."
    )
    parser.add_argument(
        "--verbosity", "-v", type=int, default=0,
        help="Conformator verbosity level (default: 0)."
    )
    parser.add_argument(
        "--output-dir", "-o", default=None,
        help="Directory for output (default: derived from input basename)."
    )
    parser.add_argument(
        "--keep-temp", action="store_true",
        help="Keep the intermediate working directory instead of removing it."
    )
    return parser.parse_args()



def run(cmd: str, description: str = "") -> None:
    """Run a shell command, exiting on failure."""
    if description:
        print(f"  {description}")
    print(f"  $ {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        print(f"ERROR: command failed (exit {ret}): {cmd}", file=sys.stderr)
        sys.exit(ret)


def main() -> None:
    args = parse_args()

    input_path = Path(args.input).resolve()
    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    basename = input_path.stem
    work_dir = Path(args.output_dir).resolve() if args.output_dir else Path.cwd() / basename

    # Set up working directory
    work_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(input_path), str(work_dir / input_path.name))
    os.chdir(str(work_dir))

    sdf_name = input_path.name

    # Activate conformator license
    run(f"{args.conformator_path} --license {args.license_key}",
        "Activating conformator license...")

    # Build conformator flags
    flags = []
    if args.keep_3d:
        flags.append("--keep3d")
    if args.hydrogens:
        flags.append("--hydrogens")

    # Run conformator
    run(
        f"{args.conformator_path} -i {sdf_name} -o confs.sdf"
        f" {' '.join(flags)} -n {args.num_conformers} -v {args.verbosity}",
        "Running conformator..."
    )

    # Split into individual conformer SDFs via babel
    run("babel confs.sdf individual_conf.sdf -m",
        "Splitting conformers with babel...")

    # Prepare directories
    os.makedirs("single_conf_sdfs", exist_ok=True)
    os.makedirs("single_conf_params", exist_ok=True)

    molecule_name_dict: dict[str, int] = {}
    lig_name_list_path = f"{basename}_lig_name_list.txt"

    # Process each individual conformer
    for dirpath, _dirnames, filenames in os.walk(os.getcwd()):
        for filename in filenames:
            if not (filename.startswith("individual_conf") and filename.endswith(".sdf")):
                continue

            # Read first line (molecule name)
            with open(filename, "r") as fh:
                molecule_name = fh.readline().strip()

            if not molecule_name:
                continue

            # Assign unique ID
            molecule_name_dict[molecule_name] = molecule_name_dict.get(molecule_name, 0) + 1
            unique_name = f"{molecule_name}_{molecule_name_dict[molecule_name]}"

            # Rewrite SDF with unique name on first line
            with open(filename, "r") as fh_in, open(f"{unique_name}.sdf", "w") as fh_out:
                first_line = True
                for line in fh_in:
                    if first_line:
                        fh_out.write(f"{unique_name}\n")
                        first_line = False
                    else:
                        fh_out.write(line)

            # Append to master confs list and ligand name list
            os.system(f"cat {unique_name}.sdf >> confs_named.sdf")
            os.system(f"echo {unique_name} >> {lig_name_list_path}")

            # Generate Rosetta params
            run(
                f"python {args.molfile_to_params} {unique_name}.sdf"
                f" -n {unique_name} --keep-names --long-names --clobber --no-pdb",
                f"Generating params for {unique_name}..."
            )

            # Move outputs
            shutil.move(f"{unique_name}.sdf", f"single_conf_sdfs/{unique_name}.sdf")
            params_file = f"{unique_name}.params"
            if os.path.exists(params_file):
                shutil.move(params_file, f"single_conf_params/{params_file}")

            # Clean up intermediate
            os.remove(filename)

    # Remove unsplit conformer file
    if os.path.exists("confs.sdf"):
        os.remove("confs.sdf")

    # Tar results and clean up
    os.chdir("..")
    output_tar = f"{basename}.tar.gz"
    run(f"tar -czvf {output_tar} {basename}",
        f"Archiving results to {output_tar}...")

    if not args.keep_temp:
        shutil.rmtree(str(work_dir))
        print(f"Cleaned up working directory: {work_dir}")

    print(f"\nDone. Output: {Path.cwd() / output_tar}")


if __name__ == "__main__":
    main()