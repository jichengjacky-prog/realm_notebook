# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7: Score round 2 placements → per-batch scores_round2.csv
#
# Run:  snakemake -s workflows/step7_score2.smk \
#           --configfile yaml/config_state3.yaml --profile profile/lsf
# ═══════════════════════════════════════════════════════════════════════════════

import os

configfile: "config.yaml"
include: "shared_config.smk"

# ── Discover batches from step 6 output ──────────────────────────────────
BATCH_IDS = discover_batch_ids("batch_{batch_id}", "rosetta_round2.done")

STEP7_DONE = os.path.join(OUTPUT_DIR, ".step7_score2.done")


# ── Rules ─────────────────────────────────────────────────────────────────

rule all:
    input: STEP7_DONE


rule score_round2:
    """Score placements from round 2 and write per-batch scores CSV."""
    input:
        done_flag   = os.path.join(TMP_ROOT, "batch_{batch_id}", "rosetta_round2.done"),
    output:
        scores_csv  = os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round2.csv"),
    resources:
        mem_mb=2000,
        cpus=1,
        queue=LSF_QUEUE_DEFAULT,
        walltime=LSF_WALLTIME_DEFAULT,
    params:
        weights_file = WEIGHTS_FILE,
        realm        = REALM_LOCATION,
        batches_dir  = BATCHES_DIR,
        python_bin   = PYTHON_BIN,
    log:
        os.path.join(TMP_ROOT, "batch_{batch_id}", "score_round2.log"),
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

scores = score_placements(os.path.join('$BATCH_DIR', 'round2'), weights_file='{params.weights_file}')
motif_ratios = {{}}
for csv_path in glob.glob(os.path.join('$BATCH_DIR', 'round2', '*', '*', 'weighted_scores.csv')):
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
print(f'Scored {{len(best)}} unique ligands (from {{len(scores)}} placements) (round 2)')
" > {log} 2>&1

        # Keep test_params (conformer params), delete heavy per-residue Rosetta output
        if [ -d "$BATCH_DIR"/round2 ]; then
            for LIG_DIR in "$BATCH_DIR"/round2/*/; do
                [ -d "$LIG_DIR" ] || continue
                find "$LIG_DIR" -mindepth 1 -maxdepth 1 ! -name 'test_params' -exec rm -rf {} + 2>/dev/null
            done
        fi
        rm -f "$BATCH_DIR"/rosetta_round2_*.log
        """


rule step7_sentinel:
    """Aggregate sentinel: depends on all score_round2 completing."""
    localrule: True
    input:
        expand(os.path.join(TMP_ROOT, "batch_{batch_id}", "scores_round2.csv"),
               batch_id=BATCH_IDS),
    output:
        STEP7_DONE,
    resources:
        mem_mb=500,
        cpus=1,
    shell:
        """
        touch {output}
        echo "Step 7 (score_round2) complete"
        """
