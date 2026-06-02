"""
Extract a target molecule by ID from SDF files inside a tar archive,
generate conformers/params, and optionally compare each conformer
against a shapedb database via NNSearch.

Usage:
    # Extract + conformers only:
    python extract_sdf_replace_missing.py --tar <input.tar.gz> --target-id <ID>
                                         --license-key <key>

    # Extract + conformers + shapedb comparison:
    python extract_sdf_replace_missing.py --tar <input.tar.gz> --target-id <ID>
                                         --license-key <key> --shapedb
                                         --shapedb-chunk <00000> [--shapedb-data /path/to/enamine]
"""

import argparse
import heapq
import os
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

from rdkit import Chem

SCRIPT_DIR = Path(__file__).resolve().parent
CONFORMER_SCRIPT = SCRIPT_DIR / "make_conformers_and_params.py"
DEFAULT_ENAMINE_ROOT = "/pi/summer.thyme-umw/enamine-REAL-2.6billion"
DEFAULT_REALM_LOCATION = "/pi/summer.thyme-umw/Ji_rosetta_discovery"
SHAPEDB_SIF = "sif/shapedb_container.sif"
SHAPEDB_K = 100000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract a molecule from a tar of SDFs, generate conformers, "
                    "and optionally compare against shapedb."
    )
    # --- Extraction ---
    parser.add_argument(
        "--tar", required=True,
        help="Path to the input tar (or tar.gz) file containing SDF files."
    )
    parser.add_argument(
        "--target-id", required=True,
        help="Target molecule ID to search for (matches _Name or idnumber property)."
    )
    parser.add_argument(
        "--output-dir", default=os.getcwd(),
        help="Directory to write output files (default: current working directory)."
    )

    # --- Conformer generation ---
    parser.add_argument(
        "--license-key", default=None,
        help="Conformator license key (required unless --no-conformers)."
    )
    parser.add_argument(
        "--no-conformers", action="store_true",
        help="Skip conformer/params generation; only extract the SDF."
    )

    # --- Shapedb comparison ---
    parser.add_argument(
        "--shapedb", action="store_true",
        help="After conformer generation, compare each conformer against shapedb."
    )
    parser.add_argument(
        "--shapedb-chunk", default=None,
        help="5-digit chunk code to search against (e.g. '49750'). "
             "If not given, searches all 53085 chunks."
    )
    parser.add_argument(
        "--shapedb-data", default=DEFAULT_ENAMINE_ROOT,
        help=f"Root of enamine REAL library (default: {DEFAULT_ENAMINE_ROOT})."
    )
    parser.add_argument(
        "--shapedb-top-n", type=int, default=10,
        help="Report top N shapedb matches per conformer (default: 10)."
    )
    parser.add_argument(
        "--shapedb-sif", default=None,
        help="Path to shapedb_container.sif (default: <realm>/sif/shapedb_container.sif)."
    )
    parser.add_argument(
        "--realm-location", default=DEFAULT_REALM_LOCATION,
        help=f"Root of the rosetta_discovery project (default: {DEFAULT_REALM_LOCATION})."
    )
    return parser.parse_args()


def extract_sdf_from_tar(tar_path: str, target_id: str, output_dir: str) -> Path:
    """Extract a single molecule matching target_id from SDFs inside a tar archive."""
    output_sdf = Path(output_dir) / f"{target_id}.sdf"

    # Open tar (supports .tar, .tar.gz, .tgz)
    mode = "r:gz" if str(tar_path).endswith((".gz", ".tgz")) else "r"
    with tarfile.open(tar_path, mode) as tar:
        sdf_members = [m for m in tar.getmembers() if m.name.endswith(".sdf")]

        if not sdf_members:
            print(f"ERROR: No .sdf files found in {tar_path}", file=sys.stderr)
            sys.exit(1)

        print(f"Found {len(sdf_members)} SDF file(s) in archive. Searching for '{target_id}'...")

        with tempfile.TemporaryDirectory() as tmpdir:
            for member in sdf_members:
                tar.extract(member, path=tmpdir)
                extracted_path = Path(tmpdir) / member.name

                supplier = Chem.SDMolSupplier(str(extracted_path), removeHs=False)
                for mol in supplier:
                    if mol is None:
                        continue
                    mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
                    idnumber = mol.GetProp("idnumber") if mol.HasProp("idnumber") else None

                    if target_id in {mol_name, idnumber}:
                        writer = Chem.SDWriter(str(output_sdf))
                        writer.write(mol)
                        writer.close()
                        print(f"Extracted '{target_id}' -> {output_sdf}")
                        return output_sdf

                # Remove extracted file to free space before processing next member
                os.remove(extracted_path)

    print(f"ERROR: Could not find '{target_id}' in any SDF within {tar_path}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
#  Step 2: generate conformers + params
# ---------------------------------------------------------------------------
def run_conformer_generation(
    sdf_path: Path, license_key: str, output_dir: str
) -> Path:
    """Run make_conformers_and_params.py on the extracted SDF.
    Returns the path to the resulting .tar.gz archive.
    """
    if not CONFORMER_SCRIPT.exists():
        print(f"ERROR: conformer script not found at {CONFORMER_SCRIPT}", file=sys.stderr)
        sys.exit(1)

    cmd = (
        f"python {CONFORMER_SCRIPT}"
        f" --input {sdf_path}"
        f" --license-key {license_key}"
        f" --output-dir {output_dir}"
    )
    print(f"Running conformer generation:\n  $ {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        print(f"ERROR: conformer generation failed with exit code {ret}", file=sys.stderr)
        sys.exit(ret)

    result_tar = Path(output_dir) / f"{sdf_path.stem}.tar.gz"
    if result_tar.exists():
        print(f"Conformer output archive: {result_tar}")
        return result_tar
    else:
        print("WARNING: expected output tar not found; checking alternate locations...")
        alt = Path.cwd() / f"{sdf_path.stem}.tar.gz"
        if alt.exists():
            return alt
        print("ERROR: cannot locate conformer output archive", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
#  Step 3: extract conformer archive and collect individual conformer SDFs
# ---------------------------------------------------------------------------
def extract_conformer_archive(tar_path: Path, output_dir: str) -> Path:
    """Extract the conformer tar.gz and return path to the single_conf_sdfs directory."""
    extract_dir = Path(output_dir) / tar_path.stem  # strips .tar.gz → basename
    if extract_dir.exists():
        shutil.rmtree(str(extract_dir))
    extract_dir.mkdir(parents=True, exist_ok=True)

    print(f"Extracting conformer archive to {extract_dir} ...")
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=str(extract_dir))
    # The archive contains a single top-level dir; find single_conf_sdfs inside
    for root, dirs, _files in os.walk(str(extract_dir)):
        if "single_conf_sdfs" in dirs:
            return Path(root) / "single_conf_sdfs"
    print("ERROR: single_conf_sdfs not found in extracted archive", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
#  Step 4: run shapedb NNSearch for a single conformer against one subchunk
# ---------------------------------------------------------------------------
def run_shapedb_search(
    conformer_sdf: Path,
    chunk_str: str,
    subchunk: int,
    enamine_root: str,
    sif_path: str,
    work_dir: Path,
) -> Path | None:
    """Run shapedb NNSearch for one conformer against one subchunk db.
    Returns path to the output text file, or None on failure.
    """
    superchunk_str = str(int(chunk_str[:3]))
    archive_name = f"condensed_params_and_db_{subchunk}"
    archive_path = (
        f"{enamine_root}/{superchunk_str}/{chunk_str}/{archive_name}.tar.gz"
    )
    db_bind = f"{archive_name}/db.db"

    # Extract db.db from the subchunk archive
    if not (work_dir / db_bind).exists():
        ret = os.system(
            f"tar -xzf {archive_path} {db_bind} -C {work_dir}"
        )
        if ret != 0:
            print(f"  WARNING: failed to extract {archive_path}", file=sys.stderr)
            return None

    output_file = work_dir / f"{chunk_str}_{subchunk}_nn.txt"
    conformer_abs = conformer_sdf.resolve()
    db_abs = (work_dir / db_bind).resolve()

    cmd = (
        f"singularity exec"
        f" --bind {db_abs}:/input/db.db"
        f" --bind {conformer_abs}:/input/{conformer_sdf.name}"
        f" {sif_path}"
        f" /pharmit/src/build/shapedb -NNSearch"
        f" -k {SHAPEDB_K}"
        f" -ligand /input/{conformer_sdf.name}"
        f" -db /input/db.db"
        f" -print"
        f" > {output_file}"
    )
    print(f"    subchunk {subchunk}: searching...")
    ret = os.system(cmd)
    if ret != 0:
        print(f"  WARNING: shapedb failed for subchunk {subchunk}", file=sys.stderr)
        return None
    return output_file


# ---------------------------------------------------------------------------
#  Step 5: parse shapedb results with a min-heap, tracking subchunk origin
# ---------------------------------------------------------------------------
def parse_shapedb_results(
    result_files: list[tuple[Path, str, int]], top_n: int
) -> list[tuple[float, str, str, int]]:
    """Parse shapedb output files, return top_n entries.
    Each entry: (score, name, chunk_str, subchunk).
    Lower score = better (shapedb distance).
    """
    heap: list[tuple[float, str, str, int]] = []  # (-score, name, chunk, subchunk)

    for fpath, chunk_str, subchunk in result_files:
        if fpath is None or not fpath.exists():
            continue
        with open(fpath, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                name, score_str = parts[0], parts[1]
                try:
                    score = float(score_str)
                except ValueError:
                    continue
                entry = (-score, name, chunk_str, subchunk)
                if len(heap) < top_n:
                    heapq.heappush(heap, entry)
                elif -score > heap[0][0]:
                    heapq.heapreplace(heap, entry)

    return [(-s, n, c, sc) for s, n, c, sc in sorted(heap, key=lambda x: x[0], reverse=True)]


# ---------------------------------------------------------------------------
#  Extract shorthand params for a hit from the enamine library
# ---------------------------------------------------------------------------
def extract_shorthand_params(
    hit_name: str,
    chunk_str: str,
    subchunk: int,
    enamine_root: str,
    output_dir: str,
) -> Path | None:
    """Extract a single molecule's shorthand_params.txt from the enamine archive.
    Returns path to the extracted file, or None on failure.
    """
    superchunk_str = str(int(chunk_str[:3]))
    archive_name = f"condensed_params_and_db_{subchunk}"
    archive_path = f"{enamine_root}/{superchunk_str}/{chunk_str}/{archive_name}.tar.gz"

    # hit_name looks like "PV-006534976486_1" (ligand_conformer)
    # The shorthand params file is named "<ligand>_shorthand_params.txt"
    base_ligand = "_".join(hit_name.split("_")[:-1])  # strip trailing _N
    member_path = f"{archive_name}/single_conf_params/{base_ligand}_shorthand_params.txt"

    output_path = Path(output_dir) / f"{base_ligand}_shorthand_params.txt"

    cmd = f"tar -xzf {archive_path} -O '{member_path}' > {output_path} 2>/dev/null"
    ret = os.system(cmd)
    if ret != 0 or not output_path.exists() or output_path.stat().st_size == 0:
        # Try alternative: maybe the name has no conformer suffix
        base_ligand2 = hit_name.split("_")[0]
        member_path2 = f"{archive_name}/single_conf_params/{base_ligand2}_shorthand_params.txt"
        output_path2 = Path(output_dir) / f"{base_ligand2}_shorthand_params.txt"
        cmd2 = f"tar -xzf {archive_path} -O '{member_path2}' > {output_path2} 2>/dev/null"
        ret2 = os.system(cmd2)
        if ret2 != 0 or not output_path2.exists() or output_path2.stat().st_size == 0:
            print(f"  WARNING: could not extract shorthand params for {hit_name}", file=sys.stderr)
            return None
        return output_path2

    return output_path


# ---------------------------------------------------------------------------
#  shapedb comparison driver
# ---------------------------------------------------------------------------
def compare_conformers_against_shapedb(
    conformer_sdfs_dir: Path,
    chunk_str: str | None,
    enamine_root: str,
    sif_path: str,
    top_n: int,
    output_dir: str,
) -> None:
    """Iterate all conformer SDFs, run shapedb NNSearch, report best matches
    and extract shorthand params for the best hit of each conformer."""
    conformer_sdfs = sorted(conformer_sdfs_dir.glob("*.sdf"))
    if not conformer_sdfs:
        print("ERROR: no conformer SDFs found to compare", file=sys.stderr)
        sys.exit(1)

    print(f"\nComparing {len(conformer_sdfs)} conformer(s) against shapedb...\n")

    chunks_to_search = (
        [chunk_str] if chunk_str else [str(i).zfill(5) for i in range(53085)]
    )

    all_results: dict[str, list[tuple[float, str, str, int]]] = {}

    for sdf_path in conformer_sdfs:
        conf_name = sdf_path.stem
        print(f"--- {conf_name} ---")

        with tempfile.TemporaryDirectory() as tmpdir:
            td = Path(tmpdir)
            result_files: list[tuple[Path, str, int]] = []

            for chunk in chunks_to_search:
                for subchunk in range(10):
                    rf = run_shapedb_search(
                        sdf_path, chunk, subchunk, enamine_root, sif_path, td
                    )
                    if rf is not None and rf.exists():
                        result_files.append((rf, chunk, subchunk))

            if not result_files:
                print(f"  WARNING: no shapedb results for {conf_name}")
                continue

            top_matches = parse_shapedb_results(result_files, top_n)
            all_results[conf_name] = top_matches

            print(f"  Top {top_n} matches for {conf_name}:")
            for rank, (score, name, _chunk, _sc) in enumerate(top_matches, 1):
                print(f"    {rank:>3}. {name:40s}  score={score:.6f}")
            print()

            # Extract shorthand params for the BEST hit
            if top_matches:
                best_score, best_name, best_chunk, best_subchunk = top_matches[0]
                sp_path = extract_shorthand_params(
                    best_name, best_chunk, best_subchunk,
                    enamine_root, output_dir,
                )
                if sp_path:
                    print(f"  Extracted shorthand params for best hit -> {sp_path}\n")

    # Write summary file
    summary_path = Path(output_dir) / "shapedb_comparison_summary.txt"
    with open(summary_path, "w") as fh:
        fh.write("conformer\ttop_hit\ttop_score\tall_hits\n")
        for conf_name, matches in all_results.items():
            top_hit = matches[0][1] if matches else "N/A"
            top_score = f"{matches[0][0]:.6f}" if matches else "N/A"
            hits_str = "; ".join(f"{n}:{s:.4f}" for s, n, _c, _sc in matches)
            fh.write(f"{conf_name}\t{top_hit}\t{top_score}\t{hits_str}\n")

    print(f"Summary written to {summary_path}")


# ---------------------------------------------------------------------------
#  main
# ---------------------------------------------------------------------------
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Resolve shapedb SIF path
    sif_path = args.shapedb_sif or os.path.join(args.realm_location, SHAPEDB_SIF)

    # Step 1: extract target SDF from tar
    output_sdf = extract_sdf_from_tar(args.tar, args.target_id, args.output_dir)

    # Step 2: generate conformers & params (unless skipped)
    if not args.no_conformers:
        if not args.license_key:
            print(
                "ERROR: --license-key is required for conformer generation "
                "(or use --no-conformers to skip).",
                file=sys.stderr,
            )
            sys.exit(1)
        conformer_tar = run_conformer_generation(
            output_sdf, args.license_key, args.output_dir
        )
    else:
        print("Skipping conformer generation (--no-conformers).")
        return

    # Step 3 & 4: optionally compare against shapedb
    if args.shapedb:
        conformer_sdfs_dir = extract_conformer_archive(conformer_tar, args.output_dir)
        compare_conformers_against_shapedb(
            conformer_sdfs_dir=conformer_sdfs_dir,
            chunk_str=args.shapedb_chunk,
            enamine_root=args.shapedb_data,
            sif_path=sif_path,
            top_n=args.shapedb_top_n,
            output_dir=args.output_dir,
        )
    else:
        print("Skipping shapedb comparison (use --shapedb to enable).")


if __name__ == "__main__":
    main()