#!/usr/bin/env python3
"""
Clean PDB files for Rosetta by removing extraneous information,
converting modified residue names, and renumbering.

Outputs both the cleaned PDB and a FASTA file of the cleaned sequence.

Original: Phil Bradley, Rhiju Das, Michael Tyka, TJ Brunette, James Thompson (Baker Lab)
Edits: Steven Combs, Sam Deluca, Jordan Willis, Rocco Moretti (Meiler Lab)
Python 3 rewrite by Ji Cheng.
"""

import os
import sys
import subprocess
import gzip
from argparse import ArgumentParser
from pathlib import Path

# ── Canonical amino acid mapping ──────────────────────────────────────────
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Modified residue → canonical mapping
MODRES = {
    'MSE': 'MET', 'SEP': 'SER', 'TPO': 'THR', 'PTR': 'TYR',
    'HYP': 'PRO', 'HSE': 'SER', 'CSE': 'CYS', 'FME': 'MET',
    'MLY': 'LYS', 'MLZ': 'LYS', 'M3L': 'LYS', 'ALY': 'LYS',
    'KCX': 'LYS', 'LLP': 'LYS', 'CXM': 'MET', 'OMT': 'MET',
    'CSO': 'CYS', 'CSD': 'CYS', 'OCS': 'CYS', 'CME': 'CYS',
    'CSS': 'CYS', 'CSX': 'CYS', 'CAS': 'CYS', 'CAF': 'CYS',
    'CCS': 'CYS', 'CMT': 'CYS', 'CYG': 'CYS', 'SMC': 'CYS',
    'SNC': 'CYS', 'CSW': 'CYS', 'OCY': 'CYS', 'BUC': 'CYS',
    'PEC': 'CYS', 'C5C': 'CYS', 'C6C': 'CYS', 'EFC': 'CYS',
    'PR3': 'CYS', 'BCS': 'CYS', 'HTI': 'CYS', 'NPH': 'CYS',
    'SCH': 'CYS', 'SVA': 'CYS', 'CY1': 'CYS', 'CY3': 'CYS',
    'CY4': 'CYS', 'CY7': 'CYS', 'CY0': 'CYS', 'CYM': 'CYS',
    'CYD': 'CYS', 'CYF': 'CYS', 'CYQ': 'CYS', 'CYR': 'CYS',
    'DCY': 'CYS', 'DYS': 'CYS', 'NYS': 'CYS', 'GT9': 'CYS',
    'S2C': 'CYS', 'SC2': 'CYS', 'SCS': 'CYS', 'BCX': 'CYS',
    'ORN': 'ALA', 'DAB': 'ALA', 'DHA': 'ALA', 'AIB': 'ALA',
    'ABA': 'ALA', 'BAL': 'ALA', 'NAL': 'ALA', 'SAR': 'GLY',
    'NLE': 'LEU', 'NVA': 'VAL', 'CIR': 'ARG', 'AGM': 'ARG',
    'HMR': 'ARG', 'HRG': 'ARG', 'ARM': 'ARG', 'ARO': 'ARG',
    'DA2': 'ARG', 'NMM': 'ARG', 'OPR': 'ARG', 'MAI': 'ARG',
    'MGG': 'ARG', 'ALG': 'ARG', 'BOR': 'ARG', 'ACL': 'ARG',
    'AAR': 'ARG', 'CLG': 'ARG', 'CLH': 'ARG', 'HIA': 'HIS',
    'HSO': 'HIS', 'OHI': 'HIS', 'DDE': 'HIS', 'PSH': 'HIS',
    'MHS': 'HIS', 'NZH': 'HIS', 'NEP': 'HIS', 'HIP': 'HIS',
    'HIC': 'HIS', '3AH': 'HIS', 'HIQ': 'HIS', 'DHI': 'HIS',
    'LYN': 'LYS', 'LYZ': 'LYS', 'DLY': 'LYS', 'M2L': 'LYS',
    'FHL': 'LYS', 'GPL': 'LYS', 'IT1': 'LYS', 'LA2': 'LYS',
    'LCK': 'LYS', 'LDH': 'LYS', 'LET': 'LYS', 'LLY': 'LYS',
    'LSO': 'LYS', 'LYM': 'LYS', 'LYP': 'LYS', 'LYR': 'LYS',
    'LYX': 'LYS', 'OBS': 'LYS', 'DLS': 'LYS', 'KST': 'LYS',
    'AKL': 'LYS', 'APK': 'LYS', 'API': 'LYS', 'BLY': 'LYS',
    'C1X': 'LYS', 'CYJ': 'LYS', 'DNL': 'LYS', 'KGC': 'LYS',
    'LCX': 'LYS', 'MCL': 'LYS', 'YCM': 'CYS',
    'DPR': 'PRO', 'H5M': 'PRO', 'HY3': 'PRO', 'HSD': 'HIS',
    'PCA': 'PRO',
    'POM': 'PRO', 'PRS': 'PRO', 'DPL': 'PRO', 'LPD': 'PRO',
    'P2Y': 'PRO', '1AB': 'PRO', '2MT': 'PRO', '4FB': 'PRO',
    'DAL': 'ALA', 'DAR': 'ARG', 'DAS': 'ASP', 'DGN': 'GLN',
    'DGL': 'GLU', 'DIL': 'ILE', 'DLE': 'LEU', 'DPN': 'PHE',
    'DSE': 'SER', 'DSN': 'SER', 'DTH': 'THR', 'DTR': 'TRP',
    'DTY': 'TYR', 'DVA': 'VAL', 'DSG': 'ASN', 'DPH': 'PHE',
    'MEA': 'PHE', 'FCL': 'PHE', 'PFF': 'PHE', 'PBF': 'PHE',
    'DAH': 'PHE', 'HPH': 'PHE', 'HPE': 'PHE', 'PHI': 'PHE',
    'PHL': 'PHE', 'PHM': 'PHE', 'PM3': 'PHE', 'PPN': 'PHE',
    'PF5': 'PHE', 'PCS': 'PHE', 'PHA': 'PHE', 'PRQ': 'PHE',
    'PSA': 'PHE', 'BIF': 'PHE', 'B2F': 'PHE', '1PA': 'PHE',
    '4PH': 'PHE', 'NFA': 'PHE', 'MTY': 'PHE',
    'FOG': 'PHE', 'FRF': 'PHE', 'HPQ': 'PHE', 'CHS': 'PHE',
    'DBY': 'TYR', 'DPQ': 'TYR', 'ESB': 'TYR', 'FLT': 'TYR',
    'FTY': 'TYR', 'IYR': 'TYR', 'MBQ': 'TYR', 'NBQ': 'TYR',
    'NIY': 'TYR', 'OTY': 'TYR', 'PTH': 'TYR', 'PTM': 'TYR',
    'PAQ': 'TYR', 'CRQ': 'TYR', '1TY': 'TYR', '2TY': 'TYR',
    '3TY': 'TYR', 'B3Y': 'TYR', 'TYS': 'TYR', 'TYY': 'TYR',
    'BTR': 'TRP', 'FTR': 'TRP', 'HTR': 'TRP', 'PAT': 'TRP',
    '1TQ': 'TRP', '4DP': 'TRP', '4FW': 'TRP', '4HT': 'TRP',
    '4IN': 'TRP', '6CW': 'TRP', '23S': 'TRP', '32S': 'TRP',
    '32T': 'TRP', 'DTR': 'TRP', 'LTR': 'TRP', 'TRO': 'TRP',
    'DLE': 'LEU', 'DNE': 'LEU', 'DNM': 'LEU',
    'FLE': 'LEU', 'LED': 'LEU', 'LEF': 'LEU', 'LNT': 'LEU',
    'MHL': 'LEU', 'MLE': 'LEU', 'MLL': 'LEU', 'MNL': 'LEU',
    'NLO': 'LEU', 'PPH': 'LEU', 'DCL': 'LEU', 'HLU': 'LEU',
    'BUG': 'LEU', 'CLE': 'LEU', 'BLE': 'LEU', '2ML': 'LEU',
    'DNG': 'LEU', 'NLN': 'LEU', 'DIL': 'ILE', 'IIL': 'ILE',
    'IML': 'ILE', 'B2I': 'ILE', 'ILX': 'ILE',
    'B2V': 'VAL', 'DIV': 'VAL', 'MNV': 'VAL', 'MVA': 'VAL',
    'OAS': 'SER', 'PG1': 'SER', 'PYR': 'SER',
    'S1H': 'SER', 'SAC': 'SER', 'SBD': 'SER', 'SBG': 'SER',
    'SBL': 'SER', 'DHL': 'SER', 'FGP': 'SER', 'GVL': 'SER',
    'MC1': 'SER', 'MIS': 'SER', 'N10': 'SER', 'NC1': 'SER',
    'BG1': 'SER', 'B3S': 'SER', 'ALO': 'THR', 'BMT': 'THR',
    'CTH': 'THR', 'OLT': 'THR', 'CRO': 'THR', 'AEI': 'THR',
    'B3D': 'ASP', 'BHD': 'ASP', 'BFD': 'ASP', 'DMK': 'ASP',
    'IAS': 'ASP', 'OHS': 'ASP', 'OXX': 'ASP', 'PHD': 'ASP',
    'AB7': 'GLU', 'B3E': 'GLU', 'CGU': 'GLU', 'GMA': 'GLU',
    'ILG': 'GLU', 'LME': 'GLU', 'MEG': 'GLU', '3MD': 'ASP',
    'ACB': 'ASP', 'ASA': 'ASP', 'ASB': 'ASP', 'ASI': 'ASP',
    'ASL': 'ASP', 'A0A': 'ASP', 'AHB': 'ASN', 'AFA': 'ASN',
    'DMH': 'ASN', 'MEN': 'ASN', 'B3X': 'ASN', 'GHG': 'GLN',
    'GLH': 'GLN', 'MGN': 'GLN', 'AME': 'MET', 'ESC': 'MET',
    'FOR': 'MET', 'MHO': 'MET', 'MME': 'MET', 'MSO': 'MET',
}


def download_pdb(pdb_id, dest_dir='.'):
    """Download a PDB file from RCSB."""
    url = 'http://www.rcsb.org/pdb/files/%s.pdb.gz' % pdb_id.upper()
    dest = os.path.join(os.path.abspath(dest_dir), '%s.pdb.gz' % pdb_id)
    cmd = 'wget --quiet %s -O %s' % (url, dest)
    print(cmd)
    subprocess.run(cmd, shell=True)
    if os.path.exists(dest):
        return dest
    print("Error: didn't download file!", file=sys.stderr)
    return None


def find_pdb(name):
    """Locate a PDB file locally, trying various extensions."""
    for suffix in ('', '.pdb', '.pdb.gz', '.pdb1.gz'):
        path = name + suffix
        if os.path.exists(path):
            return path
        path_upper = name.upper() + suffix
        if os.path.exists(path_upper):
            return path_upper
    return None


def read_pdb_lines(path):
    """Read lines from a PDB file (gzipped or plain)."""
    if path.endswith('.gz'):
        with gzip.open(path, 'rt') as fh:
            return fh.readlines()
    with open(path, 'r') as fh:
        return fh.readlines()


def clean_pdb(pdb_input, chain_id, all_chains=False, remove_chain=False,
              keep_zero_occ=False):
    """Clean a PDB and return (pdb_text, fasta_dict, stem, stats, nres)."""
    # ── Load PDB ──
    filename = find_pdb(pdb_input)
    if filename is None:
        print("File for %s not found, downloading..." % pdb_input)
        filename = download_pdb(pdb_input[:4])
        if filename is None:
            return '', {}, '', {}, 0
    else:
        print("Found PDB file at %s" % filename)

    stem = os.path.basename(filename)
    for ext in ('.gz', '.pdb1', '.pdb'):
        if stem.endswith(ext):
            stem = stem[:-len(ext)]

    lines = read_pdb_lines(filename)

    use_chains = all_chains
    if chain_id in ('_', ' '):
        chain_id = ' '

    # ── Process ──
    pdb_lines = []
    fasta = {}
    residue_buf = []
    residue_aa = ''
    count = 1
    old_resnum = '   '
    stats = {'altpos': False, 'insres': False, 'modres': False, 'misdns': False}

    for line in lines:
        if line.startswith('ENDMDL'):
            break
        if len(line) < 22:
            continue
        if not use_chains and line[21] not in chain_id:
            continue
        if line[:4] != 'ATOM' and line[:6] != 'HETATM':
            continue

        resn = line[17:20]

        # Map modified residues
        if resn in MODRES:
            orig = resn
            resn = MODRES[resn]
            line = 'ATOM  ' + line[6:17] + resn + line[20:]
            if orig == 'MSE':
                if line[12:14] == 'SE':
                    line = line[:12] + ' S' + line[14:]
                if len(line) > 76 and line[76:78] == 'SE':
                    line = line[:76] + ' S' + line[78:]
            else:
                stats['modres'] = True

        if resn not in AA_3TO1:
            continue

        resnum = line[22:27]

        # New residue → flush buffer
        if resnum != old_resnum:
            if residue_buf:
                ok = _flush(residue_buf, residue_aa, count, pdb_lines, fasta)
                if ok:
                    count += 1
                else:
                    stats['misdns'] = True
            residue_buf = []
            residue_aa = AA_3TO1[resn]

        old_resnum = resnum

        if line[26] != ' ':
            stats['insres'] = True

        if line[16] != ' ':
            stats['altpos'] = True
            if line[16] == 'A':
                line = line[:16] + ' ' + line[17:]
            else:
                continue

        if remove_chain:
            line = line[:21] + ' ' + line[22:]

        if keep_zero_occ:
            line = line[:55] + ' 1.00' + line[60:]

        residue_buf.append(line)

    # Flush final residue
    if residue_buf:
        ok = _flush(residue_buf, residue_aa, count, pdb_lines, fasta)
        if not ok:
            stats['misdns'] = True

    nres = sum(len(seq) for seq in fasta.values())
    return ''.join(pdb_lines), fasta, stem, stats, nres


def _flush(residue_buf, residue_aa, count, pdb_lines, fasta):
    """Check backbone atoms and write residue to output."""
    has_ca = has_n = has_c = False
    for line in residue_buf:
        atom = line[12:16]
        occ = float(line[55:60])
        if atom == ' CA ' and occ > 0.0:
            has_ca = True
        elif atom == ' N  ' and occ > 0.0:
            has_n = True
        elif atom == ' C  ' and occ > 0.0:
            has_c = True
    if not (has_ca and has_n and has_c):
        return False

    chain = residue_buf[0][21] if residue_buf else ' '
    for line in residue_buf:
        newnum = '%4d ' % count
        pdb_lines.append(line[:22] + newnum + line[27:])

    fasta[chain] = fasta.get(chain, '') + residue_aa
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = ArgumentParser(description='Clean PDBs for Rosetta')
    parser.add_argument('pdb', help='PDB file or 4-letter code')
    parser.add_argument('chainid', help='Chain ID(s), or "ignorechain"/"nochain"')
    parser.add_argument('--nopdbout', action='store_true', help='Skip PDB output')
    parser.add_argument('--allchains', action='store_true', help='Keep all chains')
    parser.add_argument('--removechain', action='store_true', help='Strip chain IDs')
    parser.add_argument('--keepzeroocc', action='store_true', help='Keep zero-occupancy atoms')
    args = parser.parse_args()

    chain = args.chainid
    if chain == 'ignorechain':
        args.allchains = True
    elif chain == 'nochain':
        args.removechain = True
        args.allchains = True
    if not args.allchains:
        chain = chain.upper()

    pdb_text, fasta, stem, stats, nres = clean_pdb(
        args.pdb, chain,
        all_chains=args.allchains,
        remove_chain=args.removechain,
        keep_zero_occ=args.keepzeroocc,
    )
    
    # Write outputs to the same directory as the input PDB
    out_dir = os.path.dirname(os.path.abspath(args.pdb)) if os.path.dirname(args.pdb) else '.'
    if not os.path.isabs(args.pdb) and not os.path.dirname(args.pdb):
        out_dir = '.'

    flags = [('ALT', stats['altpos']), ('INS', stats['insres']),
             ('MOD', stats['modres']), ('DNS', stats['misdns'])]
    flag_str = ' '.join(f if v else '---' for f, v in flags)
    status = 'OK' if nres > 0 else 'BAD'
    chain_display = chain if chain != ' ' else '_'
    print('%s %s %5d %s %s' % (stem, chain_display, nres, flag_str, status))

    if nres == 0:
        return

    if not args.nopdbout:
        out_pdb = os.path.join(out_dir, '%s_%s.pdb' % (stem, chain_display))
        with open(out_pdb, 'w') as fh:
            fh.write(pdb_text)
            fh.write('TER\n')

    if not args.allchains:
        for ch, seq in fasta.items():
            header = '>%s_%s' % (stem, ch)
            print(header); print(seq)
            fasta_out = os.path.join(out_dir, '%s_%s.fasta' % (stem, ch))
            with open(fasta_out, 'w') as fh:
                fh.write(header + '\n' + seq + '\n')
    else:
        seq = ''.join(fasta.values())
        header = '>%s_%s' % (stem, chain_display)
        print(header); print(seq)
        fasta_out = os.path.join(out_dir, '%s_%s.fasta' % (stem, chain_display))
        with open(fasta_out, 'w') as fh:
            fh.write(header + '\n' + seq + '\n')


if __name__ == '__main__':
    main()
