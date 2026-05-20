# Rosetta Ligand Discovery — HPC‑Optimized Workflow (LSF)

This workflow runs Rosetta’s ligand‑discovery protocol safely and efficiently on an LSF‑managed HPC cluster. It has been redesigned to:

- Protect the shared filesystem (PANFS)
- Use node‑local scratch space for all heavy I/O
- Avoid metadata storms and soft‑lockups
- Support large job arrays
- Provide deterministic, maintainable behavior

You will interact with **three scripts**:

1. `run_ligand_discovery_search_controller.py`
2. `lsf_scratch_wrapper.sh`
3. `run_ligand_discovery_search.py`

---

## 1. Controller Script  
### `run_ligand_discovery_search_controller.py`

This script:

- Expands an anchor list into individual jobs  
- Writes a `joblist.txt` file  
- Submits an LSF job array  
- Ensures each job reserves enough `/tmp` space  
- Performs a scratch‑size sanity check  

### Usage

```bash
python run_ligand_discovery_search_controller.py \
    <target.pdb> \
    "63,87,96,179" \
    <motifs.motifs> \
    <test_params_dir> \
    -2 150 -9 \
    [extra_args_file]
```

### Scratch Size

By default, the workflow uses **10 GB** of node‑local scratch.

To change this:

```bash
export SCRATCH_SIZE_GB=20
```

The controller will:

- Convert this to MB  
- Request it from LSF via `-R "rusage[tmp=X]"`  
- Ensure the motif file fits  
- Prevent jobs from running on nodes without enough `/tmp`  

### Scratch‑Size Sanity Check

The controller checks:

- The motif file size  
- The requested scratch size  
- Whether the scratch is large enough to hold Rosetta outputs  

Minimum recommended scratch:

```
max(5 GB, 2 × motif file size)
```

If the requested scratch is too small, the controller will **abort with an error**.

---

## 2. Wrapper Script  
### `lsf_scratch_wrapper.sh`

This script:

- Creates a loopback filesystem inside `/tmp`
- Sizes it according to `SCRATCH_SIZE_GB`
- Mounts it
- Exports `SCRATCH_DIR`
- Runs the worker script inside the scratch filesystem
- Automatically cleans up on exit

You do **not** run this directly — the controller calls it.

---

## 3. Worker Script  
### `run_ligand_discovery_search.py`

This script:

- Copies all inputs into scratch  
- Runs Rosetta inside Singularity  
- Processes PDB outputs  
- Scores placements  
- Compresses results  

### **Where Results Are Saved**

At the end of the job, the worker script copies only the **final, lightweight outputs** back to the directory where you originally ran the controller.

Specifically:

- `placements.tar.gz`  
- All CSV scoring files  

These are copied from the scratch filesystem back to:

```
the directory where you launched the controller script
```

This is done using:

```python
final_dest = Path.cwd()
```

This ensures:

- All heavy I/O stays inside node‑local scratch  
- Only final results touch the shared filesystem  
- The output location is predictable and matches the original workflow  

---

## 4. Submitting Jobs

Once your inputs are ready:

```bash
export SCRATCH_SIZE_GB=20   # optional
python run_ligand_discovery_search_controller.py ...
```

The controller will print the `bsub` command it submits.

Logs appear in:

```
logs/<jobid>_<index>.out
logs/<jobid>_<index>.err
```

---

## 5. Outputs

Each job produces:

- `placements.tar.gz`
- One or more CSV scoring files

These are copied back to the directory where you ran the controller.

All intermediate files remain inside the scratch filesystem and are deleted automatically.

---

## 6. Why This Workflow Is Safe for HPC

- All heavy I/O is isolated to node‑local scratch  
- PANFS only sees final compressed results  
- No wildcard file operations  
- No directory churn  
- No bind‑mount storms  
- LSF ensures nodes have enough `/tmp`  
- Scratch is automatically cleaned up  

This prevents:

- PANFS overload  
- Kernel soft‑lockups  
- NUMA migration stalls  
- Node instability  

---

## 7. Troubleshooting

### “ERROR: Requested scratch too small”
Increase:

```bash
export SCRATCH_SIZE_GB=20
```

### “Not enough /tmp space on nodes”
Reduce job concurrency or request fewer GB.

### “Singularity cannot bind file”
Ensure input paths exist and are readable.

---

## 8. File Overview

```
run_ligand_discovery_search_controller.py   # Submits LSF job array
lsf_scratch_wrapper.sh                      # Creates scratch FS, runs worker
run_ligand_discovery_search.py              # Performs Rosetta job inside scratch
joblist.txt                                 # Auto-generated job argument list
logs/                                        # Output/error logs per job
```

---

## 9. Notes

- The controller and wrapper ensure that scratch size and LSF resource requests always match.
- The worker script never touches PANFS except to copy final results.
- This workflow is safe to run at scale.

