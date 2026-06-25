# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: Aggregate final results + cleanup intermediates → top_ligands.csv
#
# Run:  snakemake -s workflows/step8_aggregate.smk \
#           --configfile yaml/config_state3.yaml
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"

# ── Discover batches from step 7 output ──────────────────────────────────
BATCH_IDS, = glob_wildcards(os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round2.csv"))
BATCH_IDS = sorted(BATCH_IDS, key=int)

TOP_LIGANDS_CSV = os.path.join(OUTPUT_DIR, "top_ligands.csv")
CLEANUP_DONE    = os.path.join(OUTPUT_DIR, "cleanup.done")


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input: CLEANUP_DONE


rule aggregate_results:
    """Aggregate all per-batch scores, deduplicate by base ligand, and write top-N CSV."""
    localrule: True
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round2.csv"),
               batch_id=BATCH_IDS),
    output:
        TOP_LIGANDS_CSV,
    params:
        max_ligands     = MAX_LIGANDS,
        min_motif_ratio = MIN_MOTIF_RATIO,
        tmp_root        = TMP_ROOT,
        python_bin      = PYTHON_BIN,
    resources:
        mem_mb=16000,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "aggregate_results.log"),
    shell:
        """
        set -e
        {params.python_bin} -c "
import heapq, csv, glob, os

scores_files = glob.glob(os.path.join('{params.tmp_root}', 'batch_*', 'scores_round2.csv'))

output_dir = os.path.dirname(os.path.dirname('{params.tmp_root}'))
batches_dir = os.path.join(output_dir, 'batches')
batch_ligands_global = set()
ligand_chunk_info = {{}}
for batch_dir in glob.glob(os.path.join('{params.tmp_root}', 'batch_*')):
    batch_name = os.path.basename(batch_dir)
    batch_id = int(batch_name.split('_')[-1])
    batch_file = os.path.join(batches_dir, f'{{batch_id // 1000:04d}}', f'batch_{{batch_id}}.txt')
    if os.path.exists(batch_file):
        with open(batch_file, 'r') as fh:
            for line in fh:
                parts = line.strip().split(',')
                if len(parts) >= 4:
                    batch_ligands_global.add(parts[0])
                    ligand_chunk_info[parts[0]] = (parts[2], parts[3])

def get_base_ligand_global(placement_name):
    for lig in sorted(batch_ligands_global, key=len, reverse=True):
        if lig in placement_name:
            return lig
    return '_'.join(placement_name.split('_')[:-1])

seen = {{}}
for sf in scores_files:
    if not os.path.exists(sf):
        continue
    with open(sf, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lig = row['ligand']
            score = float(row['score'])
            mr = float(row['real_motif_ratio'])
            if mr <= {params.min_motif_ratio}:
                continue
            base = get_base_ligand_global(lig)
            chunk, subchunk = ligand_chunk_info.get(base, ('', ''))
            if base not in seen or score > seen[base][0]:
                seen[base] = (score, lig, mr, chunk, subchunk)

heap = []
for base, (score, lig, mr, chunk, subchunk) in seen.items():
    entry = (-score, lig, score, mr, chunk, subchunk)
    if len(heap) < {params.max_ligands}:
        heapq.heappush(heap, entry)
    elif score > -heap[0][0]:
        heapq.heapreplace(heap, entry)

sorted_results = sorted(heap, reverse=True)

with open('{output}', 'w', newline='') as fh:
    writer = csv.writer(fh)
    writer.writerow(['rank', 'ligand_name', 'score', 'real_motif_ratio', 'chunk', 'subchunk'])
    for rank, (neg_score, lig, score, mr, chunk, subchunk) in enumerate(sorted_results, 1):
        writer.writerow([rank, lig, f'{{score:.6f}}', f'{{mr:.6f}}', chunk, subchunk])

print(f'Wrote {{len(sorted_results)}} top ligands to {output}')
" > {log} 2>&1
        """


rule cleanup_intermediates:
    """Remove heavy intermediate files while preserving .done, .csv, and _ready files."""
    localrule: True
    input:
        top_ligands = TOP_LIGANDS_CSV,
    output:
        CLEANUP_DONE,
    params:
        tmp_root = TMP_ROOT,
        python_bin = PYTHON_BIN,
    resources:
        mem_mb=500,
        cpus=1,
    log:
        os.path.join(OUTPUT_DIR, "cleanup_intermediates.log"),
    shell:
        """
        set -e
        {params.python_bin} -c "
import shutil, glob, os, csv

# Build set of base ligand names from top_ligands.csv
# Format: rank,ligand_name,score,real_motif_ratio,chunk,subchunk
# ligand_name is like PV-XXXX_58 (placement name with residue suffix)
top_ligand_bases = set()
with open('{input.top_ligands}', 'r') as fh:
    reader = csv.DictReader(fh)
    for row in reader:
        lig_full = row.get('ligand_name', '')
        if lig_full:
            # Strip residue suffix: PV-XXXX_58 -> PV-XXXX
            base = '_'.join(lig_full.split('_')[:-1]) if '_' in lig_full else lig_full
            top_ligand_bases.add(base)
print(f'Top ligand bases to preserve: {{len(top_ligand_bases)}}')

for batch_dir in glob.glob(os.path.join('{params.tmp_root}', 'batch_*')):
    for entry in os.listdir(batch_dir):
        epath = os.path.join(batch_dir, entry)
        if os.path.isdir(epath):
            if entry == 'round1':
                shutil.rmtree(epath, ignore_errors=True)
            elif entry == 'round2':
                # Keep test_params only for top-hit ligands
                for lig_entry in os.listdir(epath):
                    lig_path = os.path.join(epath, lig_entry)
                    if os.path.isdir(lig_path):
                        if lig_entry in top_ligand_bases:
                            # Keep test_params, delete other subdirs
                            for sub in os.listdir(lig_path):
                                sub_path = os.path.join(lig_path, sub)
                                if os.path.isdir(sub_path) and sub != 'test_params':
                                    shutil.rmtree(sub_path, ignore_errors=True)
                                elif os.path.isfile(sub_path):
                                    os.remove(sub_path)
                            print(f'  Preserved test_params for top ligand: {{lig_entry}}')
                        else:
                            shutil.rmtree(lig_path, ignore_errors=True)
                # Remove any loose files in round2/
                for f in os.listdir(epath):
                    fp = os.path.join(epath, f)
                    if os.path.isfile(fp):
                        os.remove(fp)
            elif entry != 'test_params':
                shutil.rmtree(epath, ignore_errors=True)
        elif os.path.isfile(epath):
            if not (entry.endswith('.done') or
                    entry.endswith('.csv') or
                    '_ready' in entry):
                os.remove(epath)
    print(f'  Cleaned heavy intermediates in {{os.path.basename(batch_dir)}}')

with open('{output}', 'w') as fh:
    fh.write('done\\n')
print('  All heavy intermediate files cleaned; top-hit conformer params preserved.')
" > {log} 2>&1
        """
