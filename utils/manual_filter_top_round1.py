#!/usr/bin/env python3
"""
manual_filter_top_round1.py — manual stand-in for workflows/step4_filter.smk
(multithreaded edition)

Aggregates the per-placement raw_scores.csv files under
    <tmp_root>/batch_*/round1/<ligand>/<anchor>/raw_scores.csv
and writes the top-N ligand list in the same format the pipeline uses
(ligand,score,real_motif_ratio).  Nothing is deleted — round1 temporary
results are left untouched.

Batches are processed in parallel (one worker thread per batch; --workers
controls the pool size) and a progress monitor prints the remaining queue
and the number of finished batches every --monitor-interval seconds.

Usage:
    python3 utils/manual_filter_top_round1.py \
        [--tmp-root output/snakemake_state4/tmp] \
        [--out output/snakemake_state4/top_ligands_round1.txt] \
        [--max-ligands 10] [--min-motif-ratio 0.25] \
        [--batches-dir output/snakemake_state4/batches] \
        [--workers 32] [--monitor-interval 60]
"""

import argparse
import csv
import glob
import heapq
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
REALM = os.path.dirname(os.path.dirname(HERE))  # utils -> realm root


def discover_batches_dir(tmp_root):
    """batches dir (batch ligand lists) lives next to tmp_root."""
    return os.path.join(os.path.dirname(os.path.abspath(tmp_root)), "batches")


def load_batch_ligands_global(batches_dir):
    """Collect all known base ligand names from batches/batch_*.txt files."""
    batch_ligands_global = set()
    if not batches_dir or not os.path.isdir(batches_dir):
        return batch_ligands_global
    for batch_file in glob.glob(os.path.join(batches_dir, "*", "batch_*.txt")):
        try:
            with open(batch_file, "r") as fh:
                for line in fh:
                    parts = line.strip().split(",")
                    if len(parts) >= 1 and parts[0]:
                        batch_ligands_global.add(parts[0])
        except OSError:
            continue
    return batch_ligands_global


def get_base_ligand_global(placement_name, batch_ligands_global, sorted_ligands=None):
    # Fast path: the round1 dir name IS the canonical ligand name for
    # virtually all placements — an exact set lookup avoids scanning the
    # ~2.7M-name substring list per row.
    if placement_name in batch_ligands_global:
        return placement_name
    # Slow path (merged/variant names): substring scan over the
    # pre-sorted (length desc) list, exactly like step3/step4 grouping.
    for lig in sorted_ligands if sorted_ligands is not None else sorted(batch_ligands_global, key=len, reverse=True):
        if lig in placement_name:
            return lig
    # Fallback: placement_name is already the canonical ligand name (the
    # path segment right after round1/), so keep it as-is.  The old
    # '_'.join(name.split('_')[:-1]) hack returned '' for names without
    # underscores and merged unrelated ligands into one base.
    return placement_name


def process_batch(batch_dir, min_motif_ratio, batch_ligands_global, sorted_ligands=None):
    """Aggregate one batch's round1 raw_scores.csv files.

    Returns (seen, n_rows, n_files) where seen maps base ligand →
    (score, lig, mr) — the best qualifying placement of that batch.

    score is the 'total' column of each placement row
    (== -ddg with the default weights {ddg: -1.0}); header-only files
    (zero placements) contribute no rows.

    Implementation notes (this is the hot path — it was 60x slower with
    glob.glob('round1/*/*/raw_scores.csv'), which expands 3 wildcard
    levels with per-entry stats on the shared filesystem):
      * os.scandir walks only the two directory levels we need
        (round1/<lig>/<res>), never descending into placements/.
      * header-only raw_scores.csv files (105 bytes, zero placements)
        are detected via stat size and never opened.
      * base-ligand matching uses an exact set lookup first (see
        get_base_ligand_global); the sorted substring scan is a rare
        fallback for merged/variant names.
    """
    seen = {}
    n_rows = 0
    n_files = 0
    r1 = os.path.join(batch_dir, "round1")
    try:
        with os.scandir(r1) as it:
            lig_dirs = [(e.name, e.path) for e in it if e.is_dir()]
    except OSError:
        return seen, n_rows, n_files
    for lig, lig_path in lig_dirs:
        try:
            with os.scandir(lig_path) as it2:
                res_dirs = [(e.name, e.path) for e in it2 if e.is_dir()]
        except OSError:
            continue
        for res, res_path in res_dirs:
            sf = os.path.join(res_path, "raw_scores.csv")
            try:
                if os.path.getsize(sf) <= 105:
                    n_files += 1
                    continue  # header-only marker, no placements
            except OSError:
                continue
            n_files += 1
            try:
                with open(sf, "r") as fh:
                    reader = csv.DictReader(fh)
                    for row in reader:
                        if not row.get("file"):
                            continue  # header-only marker, no placements
                        try:
                            score = float(row["total"])
                            mr = float(row["real_motif_ratio"])
                        except (KeyError, ValueError):
                            continue
                        if mr <= min_motif_ratio:
                            continue
                        n_rows += 1
                        base = get_base_ligand_global(lig, batch_ligands_global, sorted_ligands)
                        if base not in seen or score > seen[base][0]:
                            seen[base] = (score, lig, mr)
            except OSError as e:
                print(f"WARNING: cannot read {sf}: {e}", file=sys.stderr)
    return seen, n_rows, n_files


def keep_top_n(seen, max_ligands):
    heap = []
    for base, (score, lig, mr) in seen.items():
        entry = (-score, lig, score, mr)
        if len(heap) < max_ligands:
            heapq.heappush(heap, entry)
        elif score > -heap[0][0]:
            heapq.heapreplace(heap, entry)
    return heap


def write_top_list(heap, out_path):
    top_info = {}
    for _, lig, score, mr in heap:
        top_info[lig] = (score, mr)
    with open(out_path, "w") as fh:
        fh.write("ligand,score,real_motif_ratio\n")
        for lig in sorted(top_info.keys()):
            score, mr = top_info[lig]
            fh.write(f"{lig},{score:.6f},{mr:.6f}\n")
    return len(top_info)


def _batch_sort_key(path):
    """Sort batch_<id> dirs numerically when possible (batch_1000 < batch_2)."""
    name = os.path.basename(path)
    try:
        return (0, int(name.split("_")[-1]))
    except ValueError:
        return (1, name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tmp-root",
        default=os.path.join(REALM, "output", "snakemake_state4", "tmp"),
        help="directory holding batch_*/round1 data (default: output/snakemake_state4/tmp)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="output top-ligand CSV (default: <tmp_root>/../top_ligands_round1.txt)",
    )
    parser.add_argument("--max-ligands", type=int, default=10)
    parser.add_argument("--min-motif-ratio", type=float, default=0.25)
    parser.add_argument(
        "--batches-dir",
        default=None,
        help="batches dir with batch_*.txt ligand lists (default: <tmp_root>/../batches)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=min(32, (os.cpu_count() or 4) + 4),
        help="number of worker threads (default: min(32, cpu+4))",
    )
    parser.add_argument(
        "--monitor-interval",
        type=float,
        default=60.0,
        help="progress print interval in seconds (default: 60)",
    )
    args = parser.parse_args()

    tmp_root = os.path.abspath(args.tmp_root)
    out_path = args.out or os.path.join(os.path.dirname(tmp_root), "top_ligands_round1.txt")
    batches_dir = args.batches_dir or discover_batches_dir(tmp_root)

    if not os.path.isdir(tmp_root):
        print(f"ERROR: tmp root not found: {tmp_root}", file=sys.stderr)
        sys.exit(1)

    print(f"tmp root        : {tmp_root}")
    print(f"batches dir     : {batches_dir}")
    print(f"output          : {out_path}")
    print(f"max ligands     : {args.max_ligands}")
    print(f"min motif ratio : {args.min_motif_ratio}")
    print(f"workers         : {args.workers}")
    print(f"monitor interval: {args.monitor_interval}s")

    batch_ligands_global = load_batch_ligands_global(batches_dir)
    print(f"known base ligands from batches dir: {len(batch_ligands_global)}")

    # ── Enumerate batches ────────────────────────────────────────────────
    batch_dirs = sorted(
        (
            d
            for d in glob.glob(os.path.join(tmp_root, "batch_*"))
            if os.path.isdir(d)
        ),
        key=_batch_sort_key,
    )
    total_batches = len(batch_dirs)
    print(f"batches to scan: {total_batches}")
    if total_batches == 0:
        print("WARNING: no batch dirs found under tmp root — nothing to write")
        sys.exit(0)

    # ── Shared progress state (monitor thread only reads; main thread writes) ──
    state = {
        "finished": 0,
        "rows": 0,
        "n_files": 0,
        "n_seen": 0,
        "stop": threading.Event(),
    }
    started = time.monotonic()

    def monitor():
        while not state["stop"].is_set():
            state["stop"].wait(args.monitor_interval)
            if state["stop"].is_set():
                break
            elapsed = time.monotonic() - started
            remaining = total_batches - state["finished"]
            if state["finished"] > 0:
                rate = state["finished"] / elapsed  # batches per second
                eta_min = remaining / rate / 60.0 if rate > 0 else float("nan")
                eta_txt = f" | ETA: {eta_min:6.1f} min"
            else:
                eta_txt = ""
            print(
                f"[{elapsed/60.0:5.1f} min] queue: {remaining} pending | "
                f"finished batches: {state['finished']}/{total_batches} | "
                f"rows: {state['rows']} | unique ligands: {state['n_seen']}{eta_txt}",
                flush=True,
            )

    monitor_thread = threading.Thread(target=monitor, daemon=True)
    monitor_thread.start()

    # ── Parallel scan: one task per batch ────────────────────────────────
    seen = {}
    n_rows = 0
    n_files = 0

    def merge_result(future):
        nonlocal n_rows, n_files
        try:
            batch_seen, batch_rows, batch_files = future.result()
        except Exception as e:
            # One bad batch must not kill the whole scan: count it as
            # finished and keep going (its data is simply not aggregated).
            print(f"WARNING: batch failed and was skipped: {e}", file=sys.stderr)
            state["finished"] += 1
            return
        n_rows += batch_rows
        n_files += batch_files
        for base, entry in batch_seen.items():
            if base not in seen or entry[0] > seen[base][0]:
                seen[base] = entry
        state["finished"] += 1
        state["rows"] = n_rows
        state["n_files"] = n_files
        state["n_seen"] = len(seen)

    # Pre-sort the (huge) ligand set ONCE — sorting it per placement row
    # (as a naive get_base_ligand_global would) costs ~0.5s per row.
    sorted_ligands = tuple(sorted(batch_ligands_global, key=len, reverse=True))

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                process_batch, bd, args.min_motif_ratio,
                batch_ligands_global, sorted_ligands,
            )
            for bd in batch_dirs
        ]
        for fut in as_completed(futures):
            merge_result(fut)

    state["stop"].set()
    monitor_thread.join(timeout=2)

    elapsed = time.monotonic() - started
    print(
        f"scan finished in {elapsed/60.0:.1f} min: "
        f"{state['finished']}/{total_batches} batches, "
        f"{n_rows} placement rows from {n_files} raw_scores.csv files",
        flush=True,
    )
    print(f"unique base ligands with hits: {len(seen)}")

    if not seen:
        print("WARNING: no qualifying placements found — nothing to write")
        sys.exit(0)

    heap = keep_top_n(seen, args.max_ligands)
    n_top = write_top_list(heap, out_path)
    print(f"wrote {n_top} top ligands to {out_path}")

    print("\nTop ligands:")
    for _, lig, score, mr in sorted(heap, key=lambda e: -e[2]):
        print(f"  {lig}: score={score:.6f} real_motif_ratio={mr:.6f}")


if __name__ == "__main__":
    main()
