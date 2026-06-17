# Unified Discovery Controller

## Overview

`unified_discovery_controller.py` combines the entire discovery pipeline into a single controller that:

1. **Extracts conformer params** from the Enamine library (chunk/subchunk + condensed files)
2. **Generates 150 conformers per ligand** using Conformator and parameterizes them to Rosetta format
3. **Creates test_params directories** with both library conformers and generated conformers (~151 total per ligand)
4. **Runs Rosetta ligand discovery search** against a target PDB with all conformers
5. **Scores placements** using Rosetta metrics and filtering
6. **Maintains a top-1000 deduplicated heap** of best ligands (one conformer per ligand)

## Input Format

The script expects a shapedb results list in CSV format with 4+ columns:
```
score,ligand_conf,chunk,subchunk
-0.489485800266,PV-005633531035_4,29811,1
-0.492682218552,Z4448503755_1,01035,7
```

Where:
- `score`: Shape similarity score (not used in final ranking)
- `ligand_conf`: Ligand name + conformer number (e.g., `PV-005633531035_4`)
- `chunk`: Library chunk index (5 digits with leading zeros)
- `subchunk`: Subchunk within the tar file (0-9)

## Output

Final results: `<output-dir>/top_ligands.csv`
```
rank,ligand_name,score
1,PV-005633531035,12.345
2,Z4448503755,11.234
...
```

The output contains **one row per unique ligand** (deduplicates conformers), sorted by final Rosetta+filtering score.

## Usage

### Basic Example

```bash
python unified_discovery_controller.py \
  --shapedb-list /path/to/shapedb_results.csv \
  --target-pdb /path/to/target.pdb \
  --anchor-residues 79 \
  --motifs-file /pi/summer.thyme-umw/Ji_rosetta_discovery/motifs/FINAL_motifs_list_filtered_2_3_2023.motifs \
  --output-dir ./discovery_results \
  --num-conformers 150
```

### With Multiple Anchors & Tuning

```bash
python unified_discovery_controller.py \
  --shapedb-list results.csv \
  --target-pdb target.pdb \
  --anchor-residues "79,55,103" \
  --motifs-file motifs.txt \
  --output-dir ./results \
  --max-ligands 5000 \
  --batch-size 50 \
  --workers 8 \
  --num-conformers 150 \
  --license-key "YOUR_LICENSE_KEY" \
  --atr -2.5 --rep 100 --ddg -8
```

### Example: Different Conformer Count

```bash
python unified_discovery_controller.py \
  --shapedb-list results.csv \
  --target-pdb target.pdb \
  --anchor-residues 79 \
  --motifs-file motifs.txt \
  --output-dir ./results \
  --num-conformers 250  # Generate 250 conformers instead of default 150
```

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `--shapedb-list` | **required** | Path to shapedb CSV results |
| `--target-pdb` | **required** | Target protein PDB file |
| `--anchor-residues` | **required** | Anchor residue(s), Rosetta-indexed (e.g., `79` or `11,79,55`) |
| `--motifs-file` | **required** | Path to Rosetta motifs file |
| `--output-dir` | **required** | Directory for all results and logs |
| `--realm-location` | `/pi/summer.thyme-umw/Ji_rosetta_discovery` | Repo root |
| `--enamine-path` | `/pi/summer.thyme-umw/enamine-REAL-2.6billion` | Enamine library path |
| `--workers` | `4` | Parallel batch workers |
| `--max-ligands` | `1000` | Final heap size (top N ligands) |
| `--batch-size` | `100` | Ligands per batch job |
| `--num-conformers` | `150` | Conformers to generate per ligand (Step 1.5) |
| `--license-key` | `` | Conformator license key (optional) |
| `--atr` | `-2.0` | fa_atr Rosetta cutoff |
| `--rep` | `150.0` | fa_rep Rosetta cutoff |
| `--ddg` | `-9.0` | ddg Rosetta cutoff |

## Pipeline Flow

```
┌─ Read shapedb list ──┐
│                      ▼
│   Split into batches (default: 100/batch)
│                      │
│   For each batch (parallel with --workers):
│                      │
│                      ▼
│   Step 1: Extract library conformer params
│   from Enamine library
│   (chunk/subchunk → tar → shorthand → fixed)
│                      │
│                      ▼
│   Step 1.5: Generate 150 conformers per ligand
│   From SMILES using Conformator
│   Then parameterize to Rosetta format
│   (SMILES → conformators → molfile_to_params)
│                      │
│                      ▼
│   Step 2: Run Rosetta discovery search
│   With library + generated conformers (~151 per ligand)
│   (via singularity container)
│                      │
│                      ▼
│   Step 3: Score placements
│   (Rosetta metrics + filtering)
│                      │
│                      ▼
│   Step 4: Update thread-safe heap
│   Keep top-N unique ligands (no conf duplicates)
│                      │
└───────────────────────┘
          │
          ▼
   Write top_ligands.csv
   (final results, rank 1-N)
```

## Deduplication Logic

The heap ensures **one entry per unique ligand** (ligand name without conformer number):
- First conformer of a ligand that scores high enough enters the heap
- Later conformers of the same ligand are skipped
- This prevents a single ligand with multiple conformers from flooding the results

## Step 1.5: Conformer Generation

**What happens:**
- For each unique ligand in the batch, extracts SMILES from the library conformer
- Uses Conformator (via singularity container) to generate 150 new conformers
- Runs `molfile_to_params.py` on each conformer to create Rosetta-compatible params files
- Adds all generated conformers to the `test_params/` directory alongside library conformers

**Why it matters:**
- Provides conformational sampling beyond what's in the Enamine library
- 150 new conformers + 1 library conformer = ~151 total per ligand for discovery
- Increases the chance of finding favorable binding poses
- Conformers are parameterized in parallel within each batch

**Parameters:**
- `--num-conformers` (default: 150): Number of conformers to generate per ligand
- `--license-key` (optional): Conformator license key if required by your installation

**Container dependency:**
- Requires `sif/conformator_container.sif` in the realm directory
- The container includes Conformator 1.2.1, OBabel, and molfile_to_params.py

## Output Files

```
./discovery_results/
├── batch_00000/           # Per-batch working dirs
│   ├── test_params/
│   ├── 79/                # Per-anchor residue dir
│   │   └── args           # Rosetta args file
│   └── error_log.txt
├── batch_00001/
├── ...
├── error_log.txt          # Main error log
└── top_ligands.csv        # Final ranked results
```

## Performance Notes

- **Parallel workers**: Set `--workers` based on available CPU. Each worker runs a full discovery pipeline, which is I/O-heavy (tar extraction, Rosetta runs).
- **Batch size**: Smaller batches = more fine-grained progress tracking; larger batches = fewer Rosetta jobs submitted.
- **Max ligands**: The heap grows to this size, then only keeps better scores. Memory use is O(max_ligands).

## Error Handling

- Failures in individual batches are logged to `error_log.txt` and do not stop the pipeline
- Malformed lines in the shapedb list are skipped with a warning
- Failed conformer extractions are recorded; other ligands in the batch continue

## Integration with Existing Scripts

This controller **replaces** the need to manually run:
1. `prepare_test_params_directories_controller.py`
2. `prepare_test_params_directories_sub.py`
3. `run_ligand_discovery_search_controller.py`
4. `score_placed_ligands_with_filtering.py`

All steps are coordinated in a single invocation with a unified output.

## Example: Full Workflow

```bash
# 1. Run nnsearch_controller.py to get shapedb results (existing)
python nnsearch_controller.py \
  ligand.sdf \
  -w ./shapedb_results \
  -n 3000000 -j 500

# 2. Run unified discovery (new)
python unified_discovery_controller.py \
  --shapedb-list ./shapedb_results/combined_results/combined_best_3000000_chunks_00000_53084.txt \
  --target-pdb ./target.pdb \
  --anchor-residues 79 \
  --motifs-file ./motifs/FINAL_motifs_list_filtered_2_3_2023.motifs \
  --output-dir ./discovery_results \
  --workers 16 \
  --max-ligands 5000

# 3. View top results
head -20 ./discovery_results/top_ligands.csv
```
