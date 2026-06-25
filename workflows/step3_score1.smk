# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: Score round 1 placements → per-batch scores_round1.csv
#
# Run:  snakemake -s workflows/step3_score1.smk \
#           --configfile yaml/config_state3.yaml --profile profile/lsf
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"

# ── Discover batches from step 2 output ──────────────────────────────────
BATCH_IDS = discover_batch_ids("batch_{batch_id}", "rosetta_round1.done")

# ── Step sentinel ────────────────────────────────────────────────────────
STEP3_DONE = os.path.join(OUTPUT_DIR, ".step3_score1.done")


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input: STEP3_DONE


rule score_round1:
    """Score placements from round 1 and write per-batch scores CSV."""
    input:
        done_flag   = os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round1.done"),
    output:
        scores_csv  = os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round1.csv"),
    params:
        weights_file = WEIGHTS_FILE,
        realm        = REALM_LOCATION,
        batches_dir  = BATCHES_DIR,
        python_bin   = PYTHON_BIN,
    resources:
        mem_mb=2000,
        cpus=1,
        queue=LSF_QUEUE_DEFAULT,
        walltime=LSF_WALLTIME_DEFAULT,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "score_round1.log"),
    shell:
        """
        set -e
        BATCH_DIR=$(dirname {output.scores_csv})

        {params.python_bin} -c "
import sys, os, csv, glob
sys.path.insert(0, '{params.realm}/function/discovery')
from utils import score_placements

batch_id = int('{wildcards.batch_id}')
batch_file = os.path.join('{params.batches_dir}', f'{{batch_id // 1000:04d}}', f'batch_{{batch_id}}.txt')
batch_ligands = set()
with open(batch_file, 'r') as fh:
    for line in fh:
        parts = line.strip().split(',')
        if len(parts) >= 1:
            batch_ligands.add(parts[0])

def get_base_ligand(placement_name):
    for lig in sorted(batch_ligands, key=len, reverse=True):
        if lig in placement_name:
            return lig
    return '_'.join(placement_name.split('_')[:-1])

scores = score_placements(os.path.join('$BATCH_DIR', 'round1'), weights_file='{params.weights_file}')
motif_ratios = {{}}
for csv_path in glob.glob(os.path.join('$BATCH_DIR', 'round1', '*', '*', 'weighted_scores.csv')):
    with open(csv_path, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pdb_name = os.path.basename(row.get('file', ''))
            lig_name = pdb_name.replace('.pdb', '')
            if lig_name:
                try:
                    motif_ratios[lig_name] = float(row.get('real_motif_ratio', 0))
                except ValueError:
                    pass

best = {{}}
for lig, sc in scores.items():
    base = get_base_ligand(lig)
    if base not in best or sc > best[base][0]:
        best[base] = (sc, lig, motif_ratios.get(lig, 0.0))

with open('{output.scores_csv}', 'w', newline='') as fh:
    writer = csv.writer(fh)
    writer.writerow(['ligand', 'score', 'real_motif_ratio'])
    for base, (sc, lig, mr) in sorted(best.items(), key=lambda x: -x[1][0]):
        writer.writerow([lig, f'{{sc:.6f}}', f'{{mr:.6f}}'])
print(f'Scored {{len(best)}} unique ligands (from {{len(scores)}} placements)')
" > {log} 2>&1

        # Only clean up round1 if scores CSV was successfully written (idempotent guard)
        if [ -s {output.scores_csv} ]; then
            rm -rf "$BATCH_DIR"/round1
            rm -f "$BATCH_DIR"/rosetta_round1.log
        else
            echo "WARNING: scores_round1.csv is empty — keeping round1 data for debugging"
        fi
        """


rule step3_sentinel:
    """Aggregate sentinel: depends on all score_round1 completing."""
    localrule: True
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round1.csv"),
               batch_id=BATCH_IDS),
    output:
        STEP3_DONE,
    resources:
        mem_mb=500,
        cpus=1,
    shell:
        """
        touch {output}
        echo "Step 3 (score_round1) complete"
        """
