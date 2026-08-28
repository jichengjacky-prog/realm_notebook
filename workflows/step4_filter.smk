# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: Filter top-N hits from round 1 → top_ligands_round1.txt
#
# Reads the per-placement raw_scores.csv files directly from the round1
# residue dirs (no step-3 scores_round1.csv required) and NEVER deletes
# any round1 data — temporary results are kept for inspection/reruns.
#
# Run:  snakemake -s workflows/step4_filter.smk \
#           --configfile yaml/config_state3.yaml
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"

# ── Discover batches from step 1 output ──────────────────────────────────
BATCH_IDS, = glob_wildcards(os.path.join(TMP_ROOT, "batch_{batch_id}", ".round1_ready"))
BATCH_IDS = sorted(BATCH_IDS, key=int)


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input:
        os.path.join(OUTPUT_DIR, "top_ligands_round1.txt"),


rule filter_top_round1:
    """Aggregate round-1 raw_scores.csv files from all batches, keep top-N
    ligands via heapq.  Runs on LSF via bsub because heapq aggregation over
    all batches is memory-intensive.  Does NOT delete any round1 data."""
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", ".round1_ready"),
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
import heapq, csv, glob, os

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
    # Fallback: placement_name is already the canonical ligand name (the
    # path segment right after round1/), so keep it as-is.  The old
    # '_'.join(name.split('_')[:-1]) hack returned '' for names without
    # underscores and merged unrelated ligands into one base.
    return placement_name

# ── Collect from raw_scores.csv (per-placement rows) ─────────────────────
# raw_scores.csv lives at batch_<id>/round1/<ligand>/<anchor>/raw_scores.csv.
# Header-only files (105 bytes, zero placements) contribute no rows.
# score = 'total' column == -ddg with the default weights {{ddg: -1.0}}.
seen = {{}}
n_rows = 0
n_files = 0
for sf in glob.glob(os.path.join('{params.tmp_root}', 'batch_*', 'round1', '*', '*', 'raw_scores.csv')):
    if not os.path.exists(sf):
        continue
    n_files += 1
    # ligand name is the path segment right after round1/
    lig = sf.split('/round1/')[-1].split('/')[0]
    with open(sf, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if not row.get('file'):
                continue  # header-only marker, no placements
            try:
                score = float(row['total'])
                mr = float(row['real_motif_ratio'])
            except (KeyError, ValueError):
                continue
            if mr <= {params.min_motif_ratio}:
                continue
            n_rows += 1
            base = get_base_ligand_global(lig)
            if base not in seen or score > seen[base][0]:
                seen[base] = (score, lig, mr)

print(f'  Read {{n_rows}} placement rows from {{n_files}} raw_scores.csv files across all batches')

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

# NOTE: the round-1 best anchor per ligand is NOT recorded here; step6
# derives it from step3's scores_round1.csv (placement names embed the
# residue, "res<N>_..."), keeping step4's output to the single top list.
" > {log} 2>&1

        # NOTE: round1/ temporary data is intentionally KEPT — nothing is
        # deleted by this step. Remove it manually once results are finalized.
        if [ ! -s {output.top_list} ]; then
            echo "WARNING: top list empty — check that raw_scores.csv files exist under tmp/batch_*/round1/*/*/"
        fi
        """
