#!/usr/bin/env python3
"""
Fix a designed/truncated PDB so Rosetta can read it.

Issues fixed:
  1. Removes non-standard cap atoms: CAY, CY, OY, HY1-3 (N-term ACE caps)
     and CAT, NT, HNT, HT1-3 (C-term NME caps).  These are stripped from
     every residue that carries them — both at segment boundaries and
     wherever they appear inside the chain.

  2. Keeps the PDB as a single continuous chain (no TER records inserted).
     Rosetta will auto-detect chainbreaks where C–N distances exceed the
     peptide-bond threshold and apply its chainbreak score term.

Usage:
    python fix_designed_pdb_for_rosetta.py input.pdb output.pdb
"""

import sys

# Cap atoms that are NOT part of any standard Rosetta residue type.
# They originate from a computational-design truncation protocol and
# represent acetyl (ACE) / N-methylamide (NME) capping groups.
N_CAP_ATOMS = {'CAY', 'CY', 'OY', 'HY1', 'HY2', 'HY3'}      # N-terminal ACE cap
C_CAP_ATOMS = {'CAT', 'NT', 'HNT', 'HT1', 'HT2', 'HT3'}      # C-terminal NME cap
ALL_CAP_ATOMS = N_CAP_ATOMS | C_CAP_ATOMS


def fix_pdb(input_path, output_path):
    with open(input_path, 'r') as fh:
        lines = fh.readlines()

    cleaned = []
    removed_count = 0

    for line in lines:
        # Skip empty lines
        if not line.strip():
            continue

        record = line[:6].strip()

        if record in ('ATOM', 'HETATM'):
            atom_name = line[12:16].strip()

            # Remove every cap atom — whether it sits on a terminal residue
            # or on an internal residue inside the chain.
            if atom_name in ALL_CAP_ATOMS:
                removed_count += 1
                continue

            cleaned.append(line)

        elif record == 'TER':
            # Drop any existing TER records — we keep the chain continuous.
            # Rosetta detects chainbreaks automatically by C–N distance.
            pass

        else:
            # Pass through HEADER, REMARK, END, etc.
            cleaned.append(line)

    # Ensure final END record
    if cleaned and not cleaned[-1].strip().startswith('END'):
        cleaned.append('END\n')

    # Write output
    with open(output_path, 'w') as fh:
        for line in cleaned:
            fh.write(line)

    # Report stats
    print(f"Wrote {len(cleaned)} lines to {output_path}  "
          f"(removed {removed_count} cap-atom lines, dropped all TER records)")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.pdb output.pdb")
        sys.exit(1)
    fix_pdb(sys.argv[1], sys.argv[2])
