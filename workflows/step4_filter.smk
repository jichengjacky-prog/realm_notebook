# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: Filter top-N hits from round 1 → top_ligands_round1.txt
#
# Run:  snakemake -s workflows/step4_filter.smk \
#           --configfile yaml/config_state3.yaml
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"

# ── Discover batches from step 3 output ──────────────────────────────────
BATCH_IDS, = glob_wildcards(os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round1.csv"))
BATCH_IDS = sorted(BATCH_IDS, key=int)


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input:
        os.path.join(OUTPUT_DIR, "top_ligands_round1.txt"),


rule filter_top_round1:
    """Aggregate round-1 scores from all batches, keep top-N ligands via heapq,
    and clean up round-1 residue directories.  Runs on LSF via bsub because
    heapq aggregation over all batches is memory-intensive."""
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round1.csv"),
               batch_id=BATCH_IDS),
    output:
        top_list = os.path.join(OUTPUT_DIR, "top_ligands_round1.txt"),
    params:
        max_ligands     = MAX_LIGANDS,
        min_motif_ratio = MIN_MOTIF_RATIO,
        tmp_root        = TMP_ROOT,
        python_bin      = PYTHON_BIN,
    resources:
        mem_mb=16000,
        cpus=1,
        queue=LSF_QUEUE_DEFAULT,
        walltime=LSF_WALLTIME_DEFAULT,
    log:
        os.path.join(OUTPUT_DIR, "filter_top_round1.log"),
    shell:
        """
        set -e
        {params.python_bin} -c "
import heapq, csv, shutil, glob, os

output_dir = os.path.dirname(os.path.dirname('{params.tmp_root}'))
batches_dir = os.path.join(output_dir, 'batches')
batch_ligands_global = set()
for batch_dir in glob.glob(os.path.join('{params.tmp_root}', 'batch_*')):
    batch_name = os.path.basename(batch_dir)
    batch_id = int(batch_name.split('_')[-1])
    batch_file = os.path.join(batches_dir, f'{{batch_id // 1000:04d}}', f'batch_{{batch_id}}.txt')
    if os.path.exists(batch_file):
        with open(batch_file, 'r') as fh:
            for line in fh:
                parts = line.strip().split(',')
                if len(parts) >= 1:
                    batch_ligands_global.add(parts[0])

def get_base_ligand_global(placement_name):
    for lig in sorted(batch_ligands_global, key=len, reverse=True):
        if lig in placement_name:
            return lig
    return '_'.join(placement_name.split('_')[:-1])

seen = {{}}
for sf in glob.glob(os.path.join('{params.tmp_root}', 'batch_*', 'scores_round1.csv')):
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
            if base not in seen or score > seen[base][0]:
                seen[base] = (score, lig, mr)

heap = []
for base, (score, lig, mr) in seen.items():
    entry = (-score, lig, score, mr)
    if len(heap) < {params.max_ligands}:
        heapq.heappush(heap, entry)
    elif score > -heap[0][0]:
        heapq.heapreplace(heap, entry)

top_ligands = {{lig for _, lig, _, _ in heap}}
print(f'  Top {{len(top_ligands)}} unique ligands kept from round 1')


top_info = {{}}
for _, lig, score, mr in heap:
    top_info[lig] = (score, mr)
with open('{output.top_list}', 'w') as fh:
    fh.write('ligand,score,real_motif_ratio\\n')
    for lig in sorted(top_info.keys()):
        score, mr = top_info[lig]
        fh.write(f'{{lig}},{{score:.6f}},{{mr:.6f}}\\n')
print(f'  Wrote {{len(top_info)}} top ligands to {output.top_list}')
" > {log} 2>&1

        # round1/ dirs (weighted_scores.csv) are no longer needed once the
        # top ligand list is written. Delete them ONLY here, after step 4 —
        # nothing earlier in the pipeline may remove CSVs.
        if [ -s {output.top_list} ]; then
            rm -rf "{params.tmp_root}"/batch_*/round1
        else
            echo "WARNING: top list empty — keeping round1 data"
        fi
        """
