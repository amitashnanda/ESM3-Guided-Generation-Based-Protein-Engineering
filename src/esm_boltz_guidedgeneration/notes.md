# ESM3-Boltz Guided Generation - Usage Notes

This directory contains scripts for running ESM3 guided protein generation with Boltz-1 binding affinity scoring on HPC clusters using SLURM.

## Overview

The guided generation process:
1. Takes a protein (from PDB) and a target ligand (SMILES)
2. Masks a percentage of the protein sequence
3. Iteratively unmasks positions using ESM3 (local or Forge API)
4. Scores each candidate using Boltz-1 for ligand binding affinity
5. Selects the best candidate at each step

**ESM3 Model Options:**
- `esm3-medium-2024-08`: Uses Forge API (default) - requires token, good accuracy
- `esm3-large`: Uses Forge API - requires token, highest accuracy
- `local`: Uses esm3-open-small on local GPU(s) - free but limited accuracy

**Boltz Scoring Options:**
- `affinity`: Predicts binding affinity (log10 IC50) - for lead optimization
- `binary`: Predicts binding probability - for hit discovery
- MSA Server: Optional per-candidate MSA computation (slower but more accurate)

---

## Scripts

### `guided_generation.sh`
The main SLURM job script that runs a single guided generation job.

**Parameters (passed via environment variables):**
| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `ESM_API_TOKEN` | Yes* | - | Forge API token (*required for non-local models) |
| `SMILES` | Yes | - | SMILES string of the target ligand |
| `PDB_ID` | Yes | - | PDB identifier (e.g., `1RNT`) |
| `CHAIN_ID` | Yes | - | Chain identifier (e.g., `A`) |
| `MASKING_PCT` | No | `0.4` | Fraction of residues to mask (0.0-1.0) |
| `NUM_DECODING_STEPS` | No | `32` | Number of generation steps |
| `NUM_SAMPLES` | No | `10` | Candidates generated per step |
| `ESM_MODEL` | No | `esm3-medium-2024-08` | ESM3 model: `esm3-medium-2024-08`, `esm3-large`, or `local` |
| `CACHE_DIR` | No | `auto` | Directory for caching Boltz results |
| `SCORE_TYPE` | No | `affinity` | Score type: `affinity` or `binary` |
| `USE_MSA_SERVER` | No | `false` | Use Boltz MSA server: `true` or `false` |
| `TIME_LIMIT` | No | `15:00:00` | SLURM time limit (HH:MM:SS) |

**Note on SLURM resource defaults:** The `#SBATCH` directives in `guided_generation.sh` specify default values (`-c 8`, `--mem=16G`, `--gres=gpu:1`). When submitted via `submit_all_jobs.sh`, the `CPUS` and `MEM` columns from the config override these defaults via `sbatch` command-line flags.

**Note on log directory:** The shell script generates a `run_<PDBID>_<timestamp>` directory and passes it to `main.py` via `--run_dir` so that all output files (`.out`, `.err`, `.txt`) for the same job end up in the same directory. When submitted via `submit_all_jobs.sh`, SLURM's `--output`/`--error` are also pointed to this directory, preventing stray `slurm-*.out` files in the source directory.

---

### `submit_all_jobs.sh`
Batch submission script that reads protein-ligand pairs from a config file and submits multiple jobs. Passes per-job `CPUS` and `MEM` as `sbatch --cpus-per-task` and `--mem` overrides. Generates a `run_<PDBID>_<timestamp>` log directory per job and sets SLURM's `--output`/`--error` to point there.

**Usage:**
```bash
./submit_all_jobs.sh [config_file]
```

Default config file: `protein_ligand_pairs.conf`

---

### `protein_ligand_pairs.conf`
Configuration file defining protein-ligand pairs for batch submission.

**Format:**
```
PDB_ID,CHAIN_ID,SMILES,MASKING_PCT,NUM_DECODING_STEPS,NUM_SAMPLES,ESM_MODEL,CACHE_DIR,SCORE_TYPE,USE_MSA_SERVER,TIME_LIMIT,CPUS,MEM
```

**Example:**
```conf
# Comments start with #
# Format: PDB_ID,CHAIN_ID,SMILES,MASKING_PCT,NUM_DECODING_STEPS,NUM_SAMPLES,ESM_MODEL,CACHE_DIR,SCORE_TYPE,USE_MSA_SERVER,TIME_LIMIT,CPUS,MEM
1RNT,A,NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O,0.4,32,10,esm3-medium-2024-08,./boltz_cache/1rnt,affinity,true,12:00:00,8,16G
1DB1,A,C[C@H](CCCC(C)(C)O)...,0.3,64,15,esm3-medium-2024-08,./boltz_cache/1db1,affinity,true,20:00:00,12,32G
```

---

## Usage Examples

### Submitting Multiple Jobs (Recommended)

1. Set your ESM API token:
```bash
export ESM_API_TOKEN='your_token_here'
```

2. Edit `protein_ligand_pairs.conf` with your protein-ligand pairs

3. Run:
```bash
cd src/esm_boltz_guidedgeneration
./submit_all_jobs.sh
```

Or with a custom config file:
```bash
./submit_all_jobs.sh my_custom_pairs.conf
```

**Note:** The script will check if `ESM_API_TOKEN` is set and warn you if it's missing.

---

### Submitting a Single Job

Use `sbatch` with `--export` to pass parameters directly:

**Using default Forge API model (esm3-medium-2024-08):**
```bash
sbatch \
  --job-name="boltz_1RNT_A" \
  --time="10:00:00" \
  --export=ALL,ESM_API_TOKEN="your_token_here",SMILES="NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O",PDB_ID="1RNT",CHAIN_ID="A",MASKING_PCT="0.4",NUM_DECODING_STEPS="32",NUM_SAMPLES="10" \
  guided_generation.sh
```

**Using local model (no API token needed):**
```bash
sbatch \
  --job-name="boltz_1HII_A" \
  --time="10:00:00" \
  --export=ALL,SMILES="CC(C)[C@H](NC(C)=O)C(=O)N[C@@H](Cc1ccccc1)[C@@H](O)CN(CC2CCCCC2)NC(=O)[C@@H](NC(C)=O)C(C)C",PDB_ID="1HII",CHAIN_ID="A",MASKING_PCT="0.3",NUM_DECODING_STEPS="32",NUM_SAMPLES="10",ESM_MODEL="local" \
  guided_generation.sh
```

**With additional options (MSA server, custom cache):**
```bash
sbatch \
  --job-name="boltz_1RNT_A" \
  --time="10:00:00" \
  --export=ALL,ESM_API_TOKEN="your_token_here",SMILES="NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O",PDB_ID="1RNT",CHAIN_ID="A",MASKING_PCT="0.4",NUM_DECODING_STEPS="32",NUM_SAMPLES="10",USE_MSA_SERVER="true",CACHE_DIR="./boltz_cache/1rnt",SCORE_TYPE="affinity" \
  guided_generation.sh
```

> **Note:** `--output` and `--error` are not needed — `guided_generation.sh` automatically redirects stdout/stderr to `logs/boltz/run_<PDBID>_<timestamp>/`.

---

## Monitoring Jobs

**Check job status:**
```bash
squeue -u $USER
```

**View job output in real-time:**
```bash
# Find the latest run directory
RUN_DIR=$(ls -td /path/to/logs/boltz/run_* | head -1)
tail -f $RUN_DIR/boltz_<PDB>_<CHAIN>_<job_id>.out
```

**Cancel a job:**
```bash
scancel <job_id>
```

**Cancel all your jobs:**
```bash
scancel -u $USER
```

---

## Output Files

All output files for each run are organized in a timestamped directory structure:

### Run-specific Logs (per execution)
Each run creates a directory: `logs/boltz/run_<PDBID>_<timestamp>/`

| File | Description |
|------|-------------|
| `boltz_<PDB>_<CHAIN>_<jobid>.out` | Standard output (progress, results) |
| `boltz_<PDB>_<CHAIN>_<jobid>.err` | Standard error (warnings, errors) |
| `<PDB>_<CHAIN>_<model>_mask<X>_steps<Y>_<timestamp>.txt` | Detailed generation log |
| `gpu_usage_<jobid>.csv` | GPU utilization metrics (CSV) |
| `cpu_usage_<jobid>.csv` | CPU utilization metrics (CSV) |

### Results and Plots
`results/boltz/run_<PDBID>_<timestamp>/`

| File | Description |
|------|-------------|
| `affinity_optimization.png` | Comprehensive affinity trajectory plot with boxplots, best trajectory, and overall best |
| `timing_breakdown.png` | Step-by-step timing analysis (candidate generation, scoring, total) |
| `gpu_usage_plot.png` | GPU utilization visualization (3-panel: utilization, memory, temp/power) |
| `cpu_usage_plot.png` | CPU utilization visualization (3-panel: utilization, memory, load) |

---

## Troubleshooting

### ESM API Token Issues
If you see an error about missing `ESM_API_TOKEN`:
- Make sure you set the token when submitting: `--export=ALL,ESM_API_TOKEN="your_token"`
- Get your token from [EvolutionaryScale Forge](https://forge.evolutionaryscale.ai)
- For local model, no token is needed - set `ESM_MODEL="local"`

If you see `Illegal header value` errors:
- Your token likely has trailing whitespace. Both `guided_generation.sh` and `main.py` now auto-strip whitespace from the token, but ensure your `export` command doesn't introduce extra characters:
  ```bash
  # Correct:
  export ESM_API_TOKEN='your_token_here'
  # Wrong (trailing space):
  export ESM_API_TOKEN='your_token_here '
  ```

### ESM API Errors
If a job fails with `ESMProteinError`, the error message will be printed in the `.err` file. Common causes:
- API rate limiting
- Sequence validation issues
- Biosecurity screening (see below)

### Biosecurity Screening - Sequences of Concern

The ESM3 Forge API performs biosecurity screening and will **reject sequences flagged as potential biosecurity concerns**. This includes:
- **Viral proteins** (e.g., HIV-1 protease - PDB `1HII`)
- Toxins and other dangerous pathogens
- Select agents

**Error message:**
```
RuntimeError: ESM API error: Failure in forward_and_sample: {"status":"error","message":"We detected a potential sequence of concern. Please apply for elevated access and once access is approved you can self-disclose with potential_sequence_of_concern=True"}
```

**Solutions:**
1. **Use a different protein** that isn't flagged as a biosecurity concern
2. **Apply for elevated access** at [EvolutionaryScale Forge](https://forge.evolutionaryscale.ai) - once approved, you can add `potential_sequence_of_concern=True` to API calls

**Note:** This is an intentional safety measure by EvolutionaryScale, not a bug in the scripts.

### Job Timeout
Increase `TIME_LIMIT` in the config file or `--time` flag. Larger proteins or more decoding steps require more time.

### GPU Memory Issues
The default allocation requests 1 GPU with 80GB memory (A100). If running on smaller GPUs, reduce `NUM_SAMPLES` per step.

### Insufficient CPU/Memory (Boltz returns -inf)
If all candidates in a step return `Score: -inf` and the Boltz error log shows return code `-9` (killed by OOM), the job was allocated insufficient resources. Boltz structure prediction is memory-intensive, especially for larger proteins.

**Fix:** Increase `CPUS` and `MEM` in the config file for that protein:
```conf
# Small proteins (~100 residues): 8 CPUs, 16G is usually sufficient
1RNT,A,...,0.4,32,10,,./boltz_cache/1rnt,affinity,true,4:00:00,8,16G
# Medium proteins (~200-400 residues): 12 CPUs, 32G recommended
1DB1,A,...,0.4,32,10,,./boltz_cache/1db1,affinity,true,5:00:00,12,32G
# Large proteins (>400 residues): 16+ CPUs, 48-64G recommended
5NP8,A,...,0.4,48,10,,./boltz_cache/5np8,affinity,true,10:00:00,16,64G
```

The `CPUS` and `MEM` columns in the config override the default `#SBATCH` values (8 CPUs, 16G) set in `guided_generation.sh`.

### Monitoring Resource Usage
Each job automatically logs GPU and CPU usage metrics:
- **GPU metrics**: utilization %, memory used/total, temperature, power draw
- **CPU metrics**: utilization %, memory used/total, load average
- Logs are saved every 30 seconds in the job's `logs/boltz/run_<PDBID>_<timestamp>/` directory
- Summary statistics are printed at job completion
- Visualizations are automatically generated:
  - `gpu_usage_plot.png`: 3-panel plot showing GPU utilization, memory, temperature, and power
  - `cpu_usage_plot.png`: 3-panel plot showing CPU utilization, memory, and load average

---

## File Structure

```
esm_boltz_guidedgeneration/
├── guided_generation.sh      # SLURM job script
├── submit_all_jobs.sh        # Batch submission script
├── protein_ligand_pairs.conf # Config file for batch jobs
├── main.py                   # Main Python entry point
├── guided_generation.py      # Core generation logic
├── boltz_scoring_utils.py    # Boltz-1 scoring utilities
├── plot_resource_usage.py    # Resource usage plotting script
├── boltz_cache/              # Cached Boltz results
└── notes.md                  # This file
```
