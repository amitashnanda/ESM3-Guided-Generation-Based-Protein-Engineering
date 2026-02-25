"""
Boltz Scoring Utilities for ESM3 Guided Generation

This module provides scoring functions for protein-ligand binding affinity
prediction using Boltz-1 model during ESM3 guided generation.

Author: Adapted from Anna Su's original implementation
Date: 2025-12-06

Last updated: 2026-02-20 by Anna Su
  - Updated functions for plotting
"""

import os
import json
import subprocess
import shutil
import hashlib
from pathlib import Path
from typing import Optional, Dict, Tuple, Any, List
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from esm.sdk.api import ESMProtein
from guided_generation import GuidedDecodingScoringFunction


def get_cache_path(seq: str, smiles: str, cache_dir: str, step: int = 0, candidate_idx: int = 0) -> str:
    """Generate cache path based on sequence and SMILES hash.
    
    Returns path to a DIRECTORY that stores the full Boltz result
    (structure files, affinity JSON, confidence scores, etc.)
    
    Naming convention: step{N}_cand{M}_{seq_prefix}_{hash}
    Example: step01_cand003_AGDYVC_8e70cb71
    """
    combined = f"{seq}_{smiles}"
    h = hashlib.sha256(combined.encode()).hexdigest()
    seq_prefix = seq[:6]  # First 6 residues for quick identification
    return os.path.join(cache_dir, f"step{step:02d}_cand{candidate_idx:03d}_{seq_prefix}_{h[:8]}")


# ============================================================================
# OLD per-candidate cache path (stored only .json via joblib)
# Kept for backward compatibility with existing caches
# ============================================================================
def get_cache_path_legacy(seq: str, smiles: str, cache_dir: str) -> str:
    """Generate legacy cache path (joblib .json file) for backward compat."""
    combined = f"{seq}_{smiles}"
    h = hashlib.sha256(combined.encode()).hexdigest()
    return os.path.join(cache_dir, f"boltz_affinity_{h[:16]}.json")


def _check_cache(seq: str, smiles: str, cache_dir: str, score_type: str = "affinity") -> Tuple[Optional[float], Optional[Dict[str, Any]]]:
    """
    Check both new directory-based cache AND legacy joblib cache.
    Searches all step*/cand* dirs in cache_dir for a matching sequence hash.
    Returns (score, full_data) if found, (None, None) if not cached.
    """
    import joblib
    
    # --- NEW directory-based cache: search for matching hash in any step/cand dir ---
    combined = f"{seq}_{smiles}"
    h = hashlib.sha256(combined.encode()).hexdigest()[:8]
    # Look for any directory ending with this hash
    matching_dirs = list(Path(cache_dir).glob(f"step*_*_{h}"))
    if matching_dirs:
        cache_result_dir = str(matching_dirs[0])
    else:
        # Also check legacy naming (boltz_result_<hash>)
        legacy_hash = hashlib.sha256(combined.encode()).hexdigest()[:16]
        legacy_dir = os.path.join(cache_dir, f"boltz_result_{legacy_hash}")
        cache_result_dir = legacy_dir
    if os.path.isdir(cache_result_dir):
        # Look for affinity JSON inside predictions subdirectory
        affinity_files = list(Path(cache_result_dir).rglob("affinity_*.json"))
        if affinity_files:
            json_path = max(affinity_files, key=lambda p: p.stat().st_mtime)
            try:
                with open(json_path, "r") as jf:
                    data = json.load(jf)
                full_data = data.copy()
                if score_type == "binary":
                    affinity = float(data.get("affinity_probability_binary", None))
                else:
                    raw_affinity = float(data.get("affinity_pred_value", None))
                    affinity = -raw_affinity if raw_affinity is not None else None
                return affinity, full_data
            except Exception:
                pass
    
    # --- LEGACY joblib cache: check for old .json cache file ---
    legacy_path = os.path.abspath(get_cache_path_legacy(seq, smiles, cache_dir))
    if os.path.exists(legacy_path):
        try:
            cached_data = joblib.load(legacy_path)
            if isinstance(cached_data, tuple):
                return cached_data
            return cached_data, None
        except Exception:
            pass
    
    return None, None


def _save_to_cache(seq: str, smiles: str, cache_dir: str, boltz_result_dir: str, candidate_name: str, step: int = 0, candidate_idx: int = 0):
    """
    Save full Boltz prediction results to the cache directory.
    Copies the entire prediction folder (structure, affinity, confidence, etc.).
    
    Cache directory naming: step{N}_cand{M}_{seq_prefix}_{hash}
    
    Args:
        seq: Protein sequence (used for cache key)
        smiles: Ligand SMILES (used for cache key)
        cache_dir: Base cache directory
        boltz_result_dir: Path to the Boltz results directory (e.g., .../boltz_results_<name>)
        candidate_name: Name of the candidate (stem of the YAML file)
        step: Current generation step number
        candidate_idx: Candidate index within the step
    """
    cache_dest = os.path.abspath(get_cache_path(seq, smiles, cache_dir, step=step, candidate_idx=candidate_idx))
    
    if os.path.exists(cache_dest):
        shutil.rmtree(cache_dest, ignore_errors=True)
    
    # Copy the prediction subdirectory for this candidate
    pred_dir = os.path.join(boltz_result_dir, "predictions", candidate_name)
    if os.path.isdir(pred_dir):
        shutil.copytree(pred_dir, os.path.join(cache_dest, "predictions", candidate_name))
    
    # Also copy the MSA directory if it exists (useful for debugging)
    msa_dir = os.path.join(boltz_result_dir, "msa")
    if os.path.isdir(msa_dir):
        shutil.copytree(msa_dir, os.path.join(cache_dest, "msa"))





# ============================================================================
# BATCHED SCORING: Score all candidates in one Boltz subprocess call
# This loads the model ONCE per step instead of once per candidate (~3x faster)
# ============================================================================

def compute_boltz_affinity_batch(
    sequences: List[str],
    smiles: str,
    cache_dir: str,
    timeout_sec: int = 7200,
    cleanup_tmp: bool = True,
    verbose_boltz: bool = False,
    step: int = 1,
    log_file_path: Optional[str] = None,
    score_type: str = "affinity",
    use_msa_server: bool = False,
    num_devices: Optional[int] = None,
) -> List[Tuple[Optional[float], Optional[Dict[str, Any]]]]:
    """
    Compute protein-ligand binding affinity for MULTIPLE sequences in one Boltz call.
    
    Instead of spawning a new subprocess per candidate (which reloads the model each time),
    this creates a directory of YAML files and calls Boltz ONCE to process all of them.
    The model is loaded once and predictions run for all candidates.
    
    When use_msa_server is True, Boltz computes MSA independently for each candidate
    sequence via the --use_msa_server flag (per-candidate MSA).
    
    Args:
        sequences: List of protein sequences to score
        smiles: Ligand SMILES string (same for all candidates)
        cache_dir: Directory for caching results
        timeout_sec: Timeout for Boltz execution
        cleanup_tmp: Whether to clean up temporary files
        verbose_boltz: Whether to print verbose output
        step: Current generation step (for logging)
        log_file_path: Path to log file for errors
        score_type: "affinity" or "binary"
        use_msa_server: Whether to use MSA server for better accuracy
        num_devices: Number of GPUs to use (None = auto-detect all available)
    
    Returns:
        List of (score, full_affinity_data_dict) tuples, one per sequence.
        Score is None if prediction fails for that sequence.
    """
    results = [None] * len(sequences)
    uncached_indices = []
    uncached_sequences = []
    
    # --- Step 1: Check cache for each sequence ---
    for i, seq in enumerate(sequences):
        cached_score, cached_data = _check_cache(seq, smiles, cache_dir, score_type)
        if cached_score is not None:
            results[i] = (cached_score, cached_data)
            if verbose_boltz:
                print(f"[BOLTZ BATCH] Candidate {i+1}: cache hit (score={cached_score:.4f})")
        else:
            uncached_indices.append(i)
            uncached_sequences.append(seq)
    
    # If everything was cached, return immediately
    if not uncached_sequences:
        print(f"[BOLTZ BATCH] All {len(sequences)} candidates found in cache.")
        return results
    
    print(f"[BOLTZ BATCH] {len(sequences) - len(uncached_sequences)} cached, "
          f"{len(uncached_sequences)} need scoring.")
    
    # --- Step 2: Create a directory with YAML files for all uncached candidates ---
    # Organized as: <cache_dir>/batch_step{N}/inputs/*.yaml
    # These directories are PRESERVED (not cleaned up) for debugging/reference
    batch_dir = os.path.abspath(os.path.join(
        cache_dir, f"batch_step{step:02d}"
    ))
    yaml_dir = os.path.join(batch_dir, "inputs")
    os.makedirs(yaml_dir, exist_ok=True)
    
    candidate_names = []
    for batch_idx, (orig_idx, seq) in enumerate(zip(uncached_indices, uncached_sequences)):
        candidate_name = f"candidate_{orig_idx+1:03d}"
        candidate_names.append(candidate_name)
        
        yaml_path = os.path.join(yaml_dir, f"{candidate_name}.yaml")
        
        yaml_content = (
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: '{seq}'\n"
            "  - ligand:\n"
            "      id: B\n"
            f"      smiles: '{smiles}'\n"
            "properties:\n"
            "  - affinity:\n"
            "      binder: B\n"
        )
        with open(yaml_path, "w") as f:
            f.write(yaml_content)
    
    # --- Step 3: Run Boltz prediction on the entire directory ---
    try:
        boltz_exe = shutil.which("boltz")
        if not boltz_exe:
            raise RuntimeError("Boltz executable not found in PATH")
        
        cmd = [
            boltz_exe, "predict", yaml_dir,
            "--out_dir", batch_dir,
            "--override",
        ]
        
        # Use --use_msa_server for per-candidate MSA generation
        if use_msa_server:
            cmd.append("--use_msa_server")
        
        # Use all available GPUs via Boltz's --devices flag (DDP)
        if num_devices is None:
            num_devices = torch.cuda.device_count() if torch.cuda.is_available() else 1
        # Boltz requires devices <= num predictions
        num_devices = min(num_devices, len(uncached_sequences))
        print(f"Number of devices: {num_devices}")
        if num_devices > 1:
            cmd.extend(["--devices", str(num_devices)])
        
        # OLD: Restrict to single GPU via CUDA_VISIBLE_DEVICES
        # env = os.environ.copy()
        # if gpu_id is not None:
        #     env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        env = os.environ.copy()
        
        print(f"[BOLTZ BATCH] Running Boltz for {len(uncached_sequences)} candidates "
              f"across {num_devices} GPU(s)...")
        
        proc = subprocess.run(
            cmd,
            cwd=batch_dir,
            capture_output=True,
            text=True,
            timeout=timeout_sec,
            env=env
        )
        
        if verbose_boltz and log_file_path:
            with open(log_file_path, 'a') as f:
                f.write(f"\n--- BOLTZ BATCH Log Step {step} ---\n")
                f.write(f"STDOUT:\n{proc.stdout}\n")
                f.write(f"STDERR:\n{proc.stderr}\n")
                f.write("---------------------------------------------------\n")
        
        if proc.returncode != 0:
            # Always log stderr on failure (not just verbose)
            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write(f"[ERROR] Boltz batch failed at step {step}\n")
                    f.write(f"Return code: {proc.returncode}\n")
                    f.write(f"STDERR: {proc.stderr[:1000]}\n")
            print(f"[ERROR] Boltz batch failed with return code {proc.returncode}")
            print(f"[ERROR] STDERR: {proc.stderr[:500]}")
            # Return None for all uncached candidates
            for orig_idx in uncached_indices:
                results[orig_idx] = (None, None)
            return results
        
        # --- Step 4: Parse results for each candidate ---
        # Boltz creates: <batch_dir>/boltz_results_inputs/predictions/<candidate_name>/affinity_<candidate_name>.json
        boltz_results_dir = os.path.join(batch_dir, "boltz_results_inputs")
        
        # Check if predictions directory exists and has content
        pred_base = os.path.join(boltz_results_dir, "predictions")
        if not os.path.isdir(pred_base) or not os.listdir(pred_base):
            err_msg = (f"[ERROR] Boltz returned success but predictions directory is empty. "
                       f"This usually means MSA files are missing. "
                       f"Try adding --use_msa_server flag.")
            print(err_msg)
            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write(f"{err_msg}\n")
                    if proc.stderr:
                        f.write(f"STDERR: {proc.stderr[:1000]}\n")
            for orig_idx in uncached_indices:
                results[orig_idx] = (None, None)
            return results
        
        for batch_idx, (orig_idx, seq, cand_name) in enumerate(
            zip(uncached_indices, uncached_sequences, candidate_names)
        ):
            try:
                affinity_json = os.path.join(
                    boltz_results_dir, "predictions", cand_name, f"affinity_{cand_name}.json"
                )
                
                if os.path.exists(affinity_json):
                    with open(affinity_json, "r") as jf:
                        data = json.load(jf)
                    
                    full_data = data.copy()
                    
                    if score_type == "binary":
                        affinity = float(data.get("affinity_probability_binary", None))
                    else:
                        raw_affinity = float(data.get("affinity_pred_value", None))
                        affinity = -raw_affinity if raw_affinity is not None else None
                    
                    results[orig_idx] = (affinity, full_data)
                    
                    # Save full prediction to cache (with step/candidate info in the name)
                    # Use 1-based indexing for cache directory (cand001 instead of cand000)
                    _save_to_cache(seq, smiles, cache_dir, boltz_results_dir, cand_name, step=step, candidate_idx=orig_idx+1)
                    
                    if verbose_boltz and affinity is not None:
                        # Log message already uses 1-based indexing (orig_idx+1)
                        print(f"[BOLTZ BATCH] Candidate {orig_idx+1}: score={affinity:.4f}")
                else:
                    if log_file_path:
                        with open(log_file_path, 'a') as f:
                            f.write(f"[ERROR] No affinity JSON found for candidate {orig_idx+1} "
                                    f"(expected: {affinity_json})\n")
                    results[orig_idx] = (None, None)
                    
            except Exception as e:
                if log_file_path:
                    with open(log_file_path, 'a') as f:
                        f.write(f"[ERROR] Failed to parse result for candidate {orig_idx+1}: {e}\n")
                results[orig_idx] = (None, None)
    
    except subprocess.TimeoutExpired:
        if log_file_path:
            with open(log_file_path, 'a') as f:
                f.write(f"[ERROR] Boltz batch timeout at step {step} ({timeout_sec}s)\n")
        for orig_idx in uncached_indices:
            if results[orig_idx] is None:
                results[orig_idx] = (None, None)
    
    except Exception as e:
        if log_file_path:
            with open(log_file_path, 'a') as f:
                f.write(f"[ERROR] Exception during Boltz batch at step {step}: {e}\n")
        if verbose_boltz:
            print(f"[ERROR] Boltz batch prediction failed: {e}")
        for orig_idx in uncached_indices:
            if results[orig_idx] is None:
                results[orig_idx] = (None, None)
    
    finally:
        # OLD: Clean up the batch directory
        # Now we KEEP batch directories for debugging/reference
        # if cleanup_tmp:
        #     shutil.rmtree(batch_dir, ignore_errors=True)
        print(f"[BOLTZ BATCH] Results preserved at: {batch_dir}")
    
    return results

class BoltzScorer(GuidedDecodingScoringFunction):
    """
    Scoring function for ESM3 guided generation using Boltz-1 binding affinity.
    
    Supports two modes:
    - score_batch(): Batched scoring — all candidates scored in one Boltz call (~3x faster)
    - __call__(): Single candidate scoring (falls back to batch of 1)
    
    Higher affinity values indicate better binding.
    """
    
    def __init__(
        self,
        smiles: str,
        cache_dir: str,
        timeout_sec: int = 3600,
        cleanup_tmp: bool = True,
        verbose_boltz: bool = False,
        score_type: str = "affinity",
        use_msa_server: bool = False,
    ):
        """
        Initialize the Boltz scoring function.
        
        Args:
            smiles: SMILES string of the ligand.
            cache_dir: Directory to cache Boltz results.
            timeout_sec: Timeout for Boltz process in seconds.
            cleanup_tmp: Whether to delete temporary files after running.
            verbose_boltz: Whether to print Boltz stdout/stderr.
            score_type: "affinity" or "binary".
            use_msa_server: Whether to use Boltz MSA server (per-candidate).
        """
        self.smiles = smiles
        self.cache_dir = os.path.abspath(cache_dir)
        self.timeout_sec = timeout_sec
        self.cleanup_tmp = cleanup_tmp
        self.verbose_boltz = verbose_boltz
        self.score_type = score_type
        self.use_msa_server = use_msa_server
        
        # Define invalid amino acids (same as FoldX implementation)
        self.invalid_amino_acids = {'B', 'J', 'O', 'U', 'X', 'Z'}
    
    def score_batch(
        self,
        proteins: List[ESMProtein],
        step: int,
        log_file_path: str,
        num_devices: Optional[int] = None,
    ) -> List[Tuple[float, Optional[Dict[str, Any]]]]:
        """
        Score a batch of proteins in ONE Boltz subprocess call.
        
        This is ~3x faster than scoring individually because the Boltz model
        is loaded once for all candidates instead of once per candidate.
        
        Args:
            proteins: List of ESMProtein objects to score
            step: Current generation step
            log_file_path: Path to log file
            num_devices: Number of GPUs to use (None = auto-detect all available)
        
        Returns:
            List of (score, affinity_details_dict) tuples. Score is -inf for invalid sequences.
        """
        # Separate valid and invalid sequences
        valid_indices = []
        valid_sequences = []
        results = [(float("-inf"), None)] * len(proteins)
        
        for i, protein in enumerate(proteins):
            seq = protein.sequence
            if any(char in self.invalid_amino_acids for char in seq):
                print(f"[SCORE BATCH] Candidate {i+1}: invalid sequence (non-standard amino acid). Score: -inf")
                # results[i] already set to (-inf, None)
            else:
                valid_indices.append(i)
                valid_sequences.append(seq)
        
        if not valid_sequences:
            return results
        
        batch_results = compute_boltz_affinity_batch(
            sequences=valid_sequences,
            smiles=self.smiles,
            cache_dir=self.cache_dir,
            timeout_sec=self.timeout_sec,
            cleanup_tmp=self.cleanup_tmp,
            verbose_boltz=self.verbose_boltz,
            step=step,
            log_file_path=log_file_path,
            score_type=self.score_type,
            use_msa_server=self.use_msa_server,
            num_devices=num_devices,
        )
        
        # Map batch results back to original indices
        for batch_idx, orig_idx in enumerate(valid_indices):
            score, details = batch_results[batch_idx]
            if score is None:
                results[orig_idx] = (float("-inf"), None)
            else:
                results[orig_idx] = (float(score), details)
        
        return results
    
    def __call__(self, protein: ESMProtein, step: int, log_file_path: str, gpu_id: Optional[int] = None) -> Tuple[float, Optional[Dict[str, Any]]]:
        """
        Score a single protein. Falls back to batch scoring with batch size 1.
        
        Args:
            protein: ESMProtein object to score
            step: Current generation step
            log_file_path: Path to log file
            gpu_id: GPU ID for this worker (None = use default)
        
        Returns:
            Tuple of (score, affinity_details_dict). Score is -inf for invalid sequences.
        """
        results = self.score_batch([protein], step=step, log_file_path=log_file_path, gpu_id=gpu_id)
        return results[0]
    
    # ========================================================================
    # OLD __call__ method (per-candidate subprocess scoring)
    # Commented out — replaced by batch scoring above
    # ========================================================================
    #
    # def __call__(self, protein: ESMProtein, step: int, log_file_path: str, gpu_id: Optional[int] = None) -> Tuple[float, Optional[Dict[str, Any]]]:
    #     """
    #     Score a protein based on predicted binding affinity to the target ligand.
    #     
    #     Args:
    #         protein: ESMProtein object to score
    #         step: Current generation step
    #         log_file_path: Path to log file
    #         gpu_id: GPU ID for this worker (None = use default)
    #     
    #     Returns:
    #         Tuple of (score, affinity_details_dict). Score is -inf for invalid sequences.
    #     """
    #     sequence = protein.sequence
    #     
    #     # Check for invalid amino acids
    #     if any(char in self.invalid_amino_acids for char in sequence):
    #         print(f"[SCORE] Invalid sequence with non-standard amino acid found. Score: -inf")
    #         return float("-inf"), None
    #     
    #     # Compute binding affinity using Boltz
    #     affinity, full_data = compute_boltz_affinity(
    #         seq=sequence,
    #         smiles=self.smiles,
    #         cache_dir=self.cache_dir,
    #         timeout_sec=self.timeout_sec,
    #         cleanup_tmp=self.cleanup_tmp,
    #         verbose_boltz=self.verbose_boltz,
    #         step=step,
    #         log_file_path=log_file_path,
    #         score_type=self.score_type,
    #         use_msa_server=self.use_msa_server,
    #         gpu_id=gpu_id
    #     )
    #     
    #     if affinity is None:
    #         return float("-inf"), None
    #     
    #     # Return affinity as score + full data dict for logging
    #     # For score_type="affinity": already negated in compute_boltz_affinity, so higher = better
    #     # For score_type="binary": higher probability = better
    #     return float(affinity), full_data
    # ========================================================================


def parse_log_file(filepath):
    """
    Parses the Boltz generation log file to extract step, binding affinity, and timing.
    Uses stateful parsing because 'Step' header appears once per batch of candidates.
    
    Args:
        filepath: Path to the log file
    
    Returns:
        Tuple of (results_df, timing_df) pandas DataFrames
    """
    print("Parsing log file: {}...".format(filepath))
    
    if not os.path.exists(filepath):
        raise IOError("Log file not found: {}".format(filepath))

    data = []
    timing_data = []
    current_step = 0
    
    # Regex patterns
    # Matches: --- Step 1 Denoised Results ---
    step_pattern = re.compile(r"--- Step (\d+) Denoised Results ---")
    
    # Matches: Score: -0.8078 | Raw Boltz affinity_pred_value (log10 IC50): 0.8078
    # We want the affinity value (the second number)
    affinity_pattern = re.compile(r"Score:.*\|\s*Raw Boltz affinity_pred_value .*?:\s*([-\d\.]+)")

    # Timing patterns
    timing_header_pattern = re.compile(r"Timing for Step (\d+):")
    cand_gen_pattern = re.compile(r"Candidate Generation \(GPU\):\s*([\d\.]+)\s*s")
    scoring_pattern = re.compile(r"Scoring \(GPU/Boltz\):\s*([\d\.]+)\s*s")
    total_time_pattern = re.compile(r"Total Step Time:\s*([\d\.]+)\s*s")

    # Temporary storage
    temp_timing = {}

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            # Check for Step Header
            step_match = step_pattern.search(line)
            if step_match:
                current_step = int(step_match.group(1))
                continue

            # Check for Affinity Score
            affinity_match = affinity_pattern.search(line)
            # Only record if we are inside a valid step (step > 0 usually)
            if affinity_match:
                affinity = float(affinity_match.group(1))
                # Append data (using list of dicts)
                data.append({'step': current_step, 'affinity': affinity})
                continue

            # Timing Parsing
            timing_header_match = timing_header_pattern.search(line)
            if timing_header_match:
                temp_timing = {'step': int(timing_header_match.group(1))}
                continue
            
            if temp_timing:
                cand_match = cand_gen_pattern.search(line)
                if cand_match:
                    temp_timing['candidate_gen_time'] = float(cand_match.group(1))
                
                score_time_match = scoring_pattern.search(line)
                if score_time_match:
                    temp_timing['scoring_time'] = float(score_time_match.group(1))
                
                total_match = total_time_pattern.search(line)
                if total_match:
                    temp_timing['total_time'] = float(total_match.group(1))
                    # Complete record for this step
                    timing_data.append(temp_timing)
                    temp_timing = {} # Reset
    
    if not data:
        print("DEBUG: No data found. Please check log format.")
        raise ValueError("Could not find any Boltz results in the log file.")
        
    print("Successfully parsed {} data points.".format(len(data)))
    print("Successfully parsed {} timing records.".format(len(timing_data)))
    
    import pandas as pd
    return pd.DataFrame(data), pd.DataFrame(timing_data)


def create_optimization_plot(df, save_path, title_suffix=""):
    """
    Creates a plot showing the Affinity distribution and trajectory.
    
    Args:
        df: DataFrame with columns ['step', 'affinity']
        save_path: Path to save the plot
        title_suffix: Optional suffix for the plot title
    """
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if df.empty:
        print("Warning: Empty dataframe, cannot create plot")
        return

    print("Generating optimization plot...")
    
    # Identify best candidates (lower affinity is better for log10 IC50)
    # Group by step and find min affinity
    best_in_step_idx = df.groupby('step')['affinity'].idxmin()
    best_in_step = df.loc[best_in_step_idx]
    
    best_overall_idx = df['affinity'].idxmin()
    best_overall = df.loc[best_overall_idx]
    
    plt.style.use('seaborn-v0_8-talk')
    fig, ax = plt.subplots(figsize=(16, 9))

    # Boxplot of distribution
    sns.boxplot(x='step', y='affinity', data=df, ax=ax, palette="coolwarm", fliersize=0)
    sns.stripplot(x='step', y='affinity', data=df, ax=ax, color='black', alpha=0.4, jitter=0.2, size=5)

    # Optimization Line
    sns.lineplot(x='step', y='affinity', data=best_in_step, ax=ax,
                 color='orange', marker='o', markersize=8,
                 linewidth=3, label='Best Candidate per Step')

    # Best Star
    ax.scatter(best_overall['step'], best_overall['affinity'],
               marker='*', color='red', s=400,
               edgecolor='black', zorder=10, label='Best Overall Candidate')

    title = 'Boltz Guided Generation: Binding Affinity Trajectory'
    if title_suffix:
        title += " - {}".format(title_suffix)
    ax.set_title(title, fontsize=30, pad=20)
    ax.set_xlabel('Generation Step', fontsize=20)
    ax.set_ylabel('Predicted Affinity (log10 IC50)', fontsize=20)
    ax.legend(fontsize=12)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print("[INFO] Optimization plot saved to: {}".format(save_path))
    plt.close()


def create_timing_plot(timing_df, save_path, title_suffix=""):
    """
    Creates a plot showing timing breakdown step by step.
    
    Args:
        timing_df: DataFrame with timing columns
        save_path: Path to save the plot
        title_suffix: Optional suffix for the plot title
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if timing_df.empty:
        print("Warning: Empty timing dataframe, skipping timing plot")
        return

    print("Generating timing plot...")
    fig, ax = plt.subplots(figsize=(16, 9))
    plt.style.use('seaborn-v0_8-talk')
    
    sns.lineplot(data=timing_df, x='step', y='candidate_gen_time', marker='o', label='Candidate Generation (GPU)', ax=ax)
    sns.lineplot(data=timing_df, x='step', y='scoring_time', marker='o', label='Scoring (GPU/Boltz)', ax=ax)
    sns.lineplot(data=timing_df, x='step', y='total_time', marker='o', linestyle='--', color='black', label='Total Step Time', ax=ax)
    
    title = 'Step Execution Time'
    if title_suffix:
        title += " - {}".format(title_suffix)
    ax.set_title(title, fontsize=24, pad=20)
    ax.set_xlabel('Step', fontsize=18)
    ax.set_ylabel('Time (seconds)', fontsize=18)
    ax.legend(fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print("[INFO] Timing plot saved to: {}".format(save_path))
    plt.close()


def plot_affinity_history(log_filepath: str, output_dir: str, seq_identifier: str = ""):
    """
    Parse log file and create optimization plots.
    
    This replaces the old simple line plot with a comprehensive visualization
    showing affinity distribution, trajectory, and timing breakdown.
    
    Args:
        log_filepath: Path to the generation log file
        output_dir: Directory to save plots
        seq_identifier: Identifier for the plot title (e.g., '1RNT_A')
    """
    try:
        # Parse log file
        results_df, timing_df = parse_log_file(log_filepath)
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate affinity optimization plot
        affinity_plot_path = os.path.join(output_dir, "affinity_optimization.png")
        create_optimization_plot(results_df, affinity_plot_path, title_suffix=seq_identifier)
        
        # Generate timing plot (if timing data is available)
        if not timing_df.empty:
            timing_plot_path = os.path.join(output_dir, "timing_breakdown.png")
            create_timing_plot(timing_df, timing_plot_path, title_suffix=seq_identifier)
        
        print("\n✓ All plots generated successfully!")
        
    except Exception as e:
        print("[ERROR] Failed to generate plots: {}".format(e))
        import traceback
        traceback.print_exc()
    """
    Creates and saves a box plot of binding affinity scores over generation steps.
    
    Args:
        all_scores_history: Dictionary mapping step numbers to lists of scores
        save_path: Path to save the plot
    """
    plot_data = []
    for step, scores in all_scores_history.items():
        for score in scores:
            if score is not None and score != float('-inf'):
                plot_data.append({'step': step, 'affinity': score})
    
    if not plot_data:
        print("[WARN] No valid scores to plot.")
        return

    df = pd.DataFrame(plot_data)
    
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Create box plot
    sns.boxplot(x='step', y='affinity', data=df, ax=ax, palette="viridis")
    
    # Overlay individual points
    sns.stripplot(
        x='step', 
        y='affinity', 
        data=df, 
        ax=ax, 
        color='black', 
        alpha=0.6, 
        jitter=True
    )
    
    ax.set_title('Distribution of Predicted Binding Affinity Scores per Generation Step', fontsize=16)
    ax.set_xlabel('Generation Step', fontsize=12)
    ax.set_ylabel('Predicted Binding Affinity (pK)', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\n[INFO] Plot saved to: {save_path}")
    plt.close()