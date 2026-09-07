#!/usr/bin/env python3
"""One-time repair of legacy Rosetta inputs; see clean_false_done_flags.py.

Scan state3/state4 batch round1/round2 inputs, without touching done flags or
job scripts. --apply requires a backup directory. Exact repeated params paths
are removed in order; headers, comments and distinct conformers are preserved.
"""
import argparse
import base64
import concurrent.futures
import errno
import json
import os
from pathlib import Path
import stat
import uuid
from collections import Counter

BASE_DIR = Path(__file__).resolve().parents[1]

# Read, do not change, the process umask. Avoid an unnecessary remote chmod
# when creating the replacement already preserves all original permission bits.
try:
    PROCESS_UMASK = int(next(s for s in Path('/proc/self/status').read_text().splitlines()
                            if s.startswith('Umask:')).split()[1], 8)
except (OSError, StopIteration, ValueError):
    PROCESS_UMASK = None


def deduplicate(data):
    seen, kept, removed = set(), [], 0
    for line in data.splitlines(keepends=True):
        entry = line.strip()
        if entry and not entry.startswith(b"#") and entry.endswith(b".params"):
            if entry in seen:
                removed += 1
                continue
            seen.add(entry)
        kept.append(line)
    return b"".join(kept), removed


def same_snapshot(a, b):
    return (a.st_dev, a.st_ino, a.st_size, a.st_mtime_ns) == (b.st_dev, b.st_ino, b.st_size, b.st_mtime_ns)


def truncate_duplicate_suffix(path, before, original, repaired):
    """Shorten only an exact redundant suffix; every retained byte stays intact.

    Readers can see the old list or its complete valid prefix. No temporary
    directory entry or chmod is needed. Shared hardlinks and read-only files
    retain the replacement path so aliases/access controls are respected.
    The caller must durably back up the original before calling this function.
    """
    if before.st_nlink != 1 or not original.startswith(repaired):
        return False
    try:
        fd = os.open(path, os.O_RDWR | os.O_NOFOLLOW | os.O_NONBLOCK)
    except PermissionError:
        return False
    try:
        if not same_snapshot(before, os.fstat(fd)) or os.pread(fd, len(original) + 1, 0) != original:
            raise RuntimeError(f"Input changed during repair: {path}")
        if not same_snapshot(before, path.lstat()):
            raise RuntimeError(f"Input changed during repair: {path}")
        os.ftruncate(fd, len(repaired))
        os.fsync(fd)
    finally:
        os.close(fd)
    if path.read_bytes() != repaired:
        raise RuntimeError(f"Verification failed: {path}")
    return True


def read_snapshot(path):
    try:
        source_fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    except OSError as exc:
        if exc.errno == errno.ELOOP:
            raise ValueError(f"Not a regular input file: {path}") from exc
        raise
    with os.fdopen(source_fd, 'rb') as source:
        before = os.fstat(source.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"Not a regular input file: {path}")
        original = source.read()
    repaired, removed = deduplicate(original)
    return before, original, repaired, removed


def backup_record(path, before, original):
    return json.dumps({"path": str(path.relative_to(BASE_DIR)),
                       "mode": stat.S_IMODE(before.st_mode),
                       "original_base64": base64.b64encode(original).decode()}) + "\n"


def repair_file(path, backup_root=None, journal=None):
    before, original, repaired, removed = read_snapshot(path)
    if not removed:
        return {"files": 1}
    if backup_root is not None:
        # One durable JSONL journal per batch avoids millions of backup files.
        journal.write(backup_record(path, before, original))
        journal.flush()
        os.fsync(journal.fileno())
        return apply_backed_up_snapshot(path, before, original, repaired, removed)
    return {"files": 1, "affected": 1, "duplicates": removed}


def apply_backed_up_snapshot(path, before, original, repaired, removed):
    """Apply only after the caller has durably recorded the original bytes."""
    if removed:
        if truncate_duplicate_suffix(path, before, original, repaired):
            return {"files": 1, "affected": 1, "duplicates": removed, "suffix_truncations": 1}
        mode = stat.S_IMODE(before.st_mode)
        temporary = path.parent / ('.residue_types.repair.' + uuid.uuid4().hex)
        fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, mode)
        try:
            with os.fdopen(fd, "wb") as fh:
                fh.write(repaired)
                fh.flush()
                os.fsync(fh.fileno())
                if PROCESS_UMASK is None or mode & ~PROCESS_UMASK != mode:
                    os.fchmod(fh.fileno(), mode)
            current = path.lstat()
            if not same_snapshot(before, current):
                raise RuntimeError(f"Input changed during repair: {path}")
            os.replace(temporary, path)
            if path.read_bytes() != repaired:
                raise RuntimeError(f"Verification failed: {path}")
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    return {"files": 1, "affected": 1, "duplicates": removed}


def process_batch_buffered(batch, backup_root):
    """Back up the whole batch durably before editing any list in that batch.

    One journal fsync per batch avoids a remote sync per input while preserving
    backup-before-write ordering. Every changed input is still checked against
    its snapshot immediately before mutation, fsynced, and read back.
    """
    counts, errors, pending = Counter(), [], []
    for round_name in ("round1", "round2"):
        round_dir = batch / round_name
        if not round_dir.is_dir():
            continue
        with os.scandir(round_dir) as entries:
            for entry in entries:
                if not entry.is_dir(follow_symlinks=False):
                    continue
                path = Path(entry.path) / 'test_params' / 'residue_types.txt'
                try:
                    snapshot = read_snapshot(path)
                    if snapshot[3]:
                        pending.append((path, snapshot))
                    else:
                        counts['files'] += 1
                except FileNotFoundError:
                    continue
                except Exception as exc:
                    errors.append(f'{path}: {exc}')
    if backup_root and pending:
        journal_path = backup_root / f'{batch.parent.parent.name}_{batch.name}.jsonl'
        with journal_path.open('x') as journal:
            for path, snapshot in pending:
                journal.write(backup_record(path, snapshot[0], snapshot[1]))
            journal.flush()
            os.fsync(journal.fileno())
    for path, snapshot in pending:
        try:
            if backup_root:
                counts.update(apply_backed_up_snapshot(path, *snapshot))
            else:
                counts.update(files=1, affected=1, duplicates=snapshot[3])
        except Exception as exc:
            errors.append(f'{path}: {exc}')
    return dict(counts), errors


def process_batch(batch, backup_root):
    counts, errors = Counter(), []
    journal = None
    try:
        if backup_root:
            journal = (backup_root / f"{batch.parent.parent.name}_{batch.name}.jsonl").open("x")
        for round_name in ("round1", "round2"):
            round_dir = batch / round_name
            if not round_dir.is_dir():
                continue
            with os.scandir(round_dir) as entries:
                for entry in entries:
                    if not entry.is_dir(follow_symlinks=False):
                        continue
                    path = Path(entry.path) / "test_params" / "residue_types.txt"
                    try:
                        counts.update(repair_file(path, backup_root, journal))
                    except FileNotFoundError:
                        continue
                    except Exception as exc:
                        errors.append(f"{path}: {exc}")
    finally:
        if journal:
            journal.close()
    return dict(counts), errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--backup-dir", type=Path)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--scope", choices=("reference", "not-done", "all-inputs"), default="reference",
                        help="reference selects done-flagged batches; not-done selects batches without rosetta_round1.done")
    parser.add_argument("--batch-backup", action="store_true",
                        help="durably back up each batch before applying its verified repairs")
    args = parser.parse_args()
    if args.workers < 1 or (args.apply and args.backup_dir is None):
        parser.error("positive --workers and --backup-dir with --apply required")
    backup_root = args.backup_dir.resolve() if args.apply else None
    if backup_root:
        backup_root.mkdir(parents=True, exist_ok=False)
    batches = []
    if args.scope == "reference":
        from clean_false_done_flags import collect_batch_dirs
        batches = [Path(p) for p in collect_batch_dirs()]
    else:
        for state in ("snakemake_state3", "snakemake_state4"):
            root = BASE_DIR / "output" / state / "tmp"
            if root.is_dir():
                with os.scandir(root) as entries:
                    batches.extend(Path(e.path) for e in entries
                                   if e.name.startswith("batch_") and e.is_dir(follow_symlinks=False))
        if args.scope == 'not-done':
            batches = [p for p in batches if not (p / 'rosetta_round1.done').exists()]
    print(f"{'REPAIR' if args.apply else 'DRY RUN'}: {len(batches)} batches", flush=True)
    if backup_root:
        (backup_root / 'selection.json').write_text(json.dumps([str(p) for p in batches], indent=2) + '\n')
    total, errors = Counter(), []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        worker = process_batch_buffered if args.batch_backup else process_batch
        futures = {executor.submit(worker, p, backup_root): p for p in batches}
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            try:
                counts, failures = future.result()
                total.update(counts)
                errors.extend(failures)
            except Exception as exc:
                errors.append(f"{futures[future]}: {exc}")
            if i % 25 == 0 or i == len(batches):
                print(json.dumps({"batches": i, **total, "errors": len(errors)}), flush=True)
    report = {"applied": args.apply, "scope": args.scope, "batch_backup": args.batch_backup,
              "batches": len(batches), **total, "errors": errors}
    if backup_root:
        (backup_root / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)
    return bool(errors)


if __name__ == "__main__":
    raise SystemExit(main())
