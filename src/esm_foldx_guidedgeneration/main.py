
# All packages used in this script
import os
import uuid
import shutil
import subprocess
import hashlib
import random
import time
from typing import List, Tuple, Optional
from typing import Optional, Tuple, Dict
import torch
from tqdm import tqdm
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, GenerationConfig
from .guided_generation import ESM3GuidedDecoding, GuidedDecodingScoringFunction
from .scoring_utils import FoldXScorer, parse_pdb_chain_sequence_with_mapping, foldx_repair_pdb, plot_ddg_history
from multiprocessing import Pool, cpu_count, Manager
from functools import partial
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from multiprocessing import set_start_method
import warnings
import argparse


# --- Main Settings ---
# This ensures multiprtocessing works correctly in different environments

try:
    set_start_method('spawn')
except RuntimeError:
    pass
warnings.simplefilter(action='ignore', category=FutureWarning)

# --- Path and Execution Settings ---
FOLDX_EXEC     = "/pscratch/sd/a/ananda/ESM3-Guided-Generation-Based-Protein-Engineering/foldx/foldx_20251231"
FOLDX_WORKDIR  = "/pscratch/sd/a/ananda/ESM3-Guided-Generation-Based-Protein-Engineering/foldx"

# --- Protein and Masking Settings to test for individual protein without batch job ---

# PDB_FILENAME   = "2lwy.pdb"
# CHAIN_ID       = "A"
# MASKING_PERCENTAGE = 0.4


# --- Generation Step Settings ---

# NUM_DECODING_STEPS = 2
# NUM_SAMPLES_PER_STEP = 2
# NUM_WORKERS = 2 

# NUM_DECODING_STEPS = 32
# NUM_SAMPLES_PER_STEP = 20
# NUM_WORKERS = 20 

# --- FoldX Settings ---

NUMBER_OF_RUNS = 1
TIMEOUT_SEC    = 1800
CLEANUP_TMP    = True 
CACHE_DIR      = os.path.join(FOLDX_WORKDIR, "foldx_cache")
VERBOSE_FOLDX  = False 

USE_STARTING_MUTATIONS = False
STARTING_MUTATIONS = { 10: "F", 25: "Y", 30: "W", 41: "I", 55: "R" }


def main():
    parser = argparse.ArgumentParser(description="Run ESM3 Guided Generation with FoldX.")
    parser.add_argument(
        "--pdb_filename",
        type=str,
        required=True,
        help="The name of the PDB file located in the FOLDX_WORKDIR (e.g., '2lwy.pdb')."
    )
    parser.add_argument(
        "--chain_id",
        type=str,
        required=True,
        help="The single-letter chain ID to extract and redesign (e.g., 'A')."
    )
    parser.add_argument(
        "--masking_percentage",
        type=float,
        default=0.40,
        help="Percentage of non-frozen residues to mask (e.g., 0.4 for 40%)."
    )

    parser.add_argument(
        "--num_decoding_steps", type=int, default=32, help="Number of generation steps."
    )
    parser.add_argument(
        "--num_samples_per_step", type=int, default=20, help="Number of candidates to generate per step."
    )
    parser.add_argument(
        "--num_workers", type=int, default=20, help="Number of parallel workers for FoldX scoring."
    )

    
    args = parser.parse_args()

    PDB_FILENAME = args.pdb_filename
    CHAIN_ID = args.chain_id
    MASKING_PERCENTAGE = args.masking_percentage
    NUM_DECODING_STEPS = args.num_decoding_steps
    NUM_SAMPLES_PER_STEP = args.num_samples_per_step
    NUM_WORKERS = args.num_workers


    os.environ["TOKENIZERS_PARALLELISM"] = "false"

    # --- 1. Initial Setup ---
    os.makedirs(FOLDX_WORKDIR, exist_ok=True); os.makedirs(CACHE_DIR, exist_ok=True)
    assert os.path.isfile(os.path.join(FOLDX_WORKDIR, "rotabase.txt")), "rotabase.txt missing"
    REPAIRED_PDB_PATH = os.path.join(FOLDX_WORKDIR, f"{os.path.splitext(PDB_FILENAME)[0]}_Repair.pdb")
    if not os.path.isfile(REPAIRED_PDB_PATH): REPAIRED_PDB_PATH = foldx_repair_pdb(FOLDX_EXEC, FOLDX_WORKDIR, PDB_FILENAME,TIMEOUT_SEC)
    WT_SEQ, SEQ_TO_PDB_MAP = parse_pdb_chain_sequence_with_mapping(REPAIRED_PDB_PATH, CHAIN_ID)

     # --- 2. Create the Masked Starting Sequence ---
    if USE_STARTING_MUTATIONS:
        print("[INFO] Applying starting mutations and freezing those positions.")
        mutant_seq_list = list(WT_SEQ)
        for pos, aa in STARTING_MUTATIONS.items():
            mutant_seq_list[pos - 1] = aa
        base_sequence = "".join(mutant_seq_list)
        
        frozen_indices_0based = {i - 1 for i in STARTING_MUTATIONS.keys()}
        maskable_indices_0based = list(set(range(len(base_sequence))) - frozen_indices_0based)

    else:
        print("[INFO] Using wild-type sequence as the base for masking.")
        base_sequence = WT_SEQ
        maskable_indices_0based = list(range(len(base_sequence)))

    
    num_to_mask = int(len(maskable_indices_0based) * MASKING_PERCENTAGE)
    indices_to_mask_0based = random.sample(maskable_indices_0based, num_to_mask)
    
    # refinement_template_list = list(starting_mutant_seq)
    refinement_template_list = list(base_sequence)
    
    for i in indices_to_mask_0based: 
        refinement_template_list[i] = '_'
    refinement_start_seq = "".join(refinement_template_list)

    
    # --- Setup the Log File ---
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_filename = f"generation_log_{timestamp}.txt"
    log_dir = "../../logs"
    os.makedirs(log_dir, exist_ok=True)
    log_filepath = os.path.join(log_dir, log_filename)
    print(f"[INFO] All generated sequences will be saved to: {log_filepath}")

    header_text = (
        f"{'='*65}\n"
        f"ESM3-FoldX Guided Generation Log\n"
        f"Run Timestamp: {timestamp}\n"
        f"{'='*65}\n\n"
        f"--- RUN PARAMETERS ---\n"
        f"PDB File: {PDB_FILENAME}\n"
        f"Masking Percentage: {MASKING_PERCENTAGE * 100:.1f}%\n"
        f"Decoding Steps: {NUM_DECODING_STEPS}\n"
        f"Samples per Step: {NUM_SAMPLES_PER_STEP}\n\n"
        f"--- SEQUENCE SETUP ---\n"
        f"Wild-Type Sequence Length: {len(WT_SEQ)}\n"
        f"Original Wild-Type Sequence:\n{WT_SEQ}\n\n"
        f"Masked Template for Generation ({len(indices_to_mask_0based)} masks):\n{refinement_start_seq}\n"
        f"{'='*65}\n"
)
    with open(log_filepath, 'w') as f:
        f.write(header_text)
    print(header_text)
    
    # --- 3. Initialize Model and Run Guided Generation ---
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\nLoading ESM3 model to {device}...")
    model = ESM3.from_pretrained().to(device).float()

    scorer_kwargs = {
        'foldx_exec': FOLDX_EXEC,
        'foldx_workdir': FOLDX_WORKDIR,
        'cache_dir': CACHE_DIR, 
        'number_of_runs': 1,
        'timeout_sec': 1800,
        'cleanup_tmp': True,
        'verbose_foldx': False,
}
    scoring_function = FoldXScorer(WT_SEQ, CHAIN_ID, SEQ_TO_PDB_MAP, REPAIRED_PDB_PATH,**scorer_kwargs)
    guided_decoding = ESM3GuidedDecoding(client=model, scoring_function=scoring_function)

    print(f"\n--- Running PARALLEL Guided Generation from the MASKED template ---")
    starting_protein = ESMProtein(sequence=refinement_start_seq)
    
    # Call the built-in function, which now supports parallelism and logging
    generated_protein, all_scores, best_overall_score, best_overall_step = guided_decoding.guided_generate(
        protein=starting_protein,
        num_decoding_steps=NUM_DECODING_STEPS,
        num_samples_per_step=NUM_SAMPLES_PER_STEP,
        track="sequence",
        num_workers=NUM_WORKERS,
        log_file_path=log_filepath
    )
    
    # --- 4. Display and Save Final Result ---
    if generated_protein:
        # best_last_step_score = max(s for s in all_scores.get(NUM_DECODING_STEPS, []) if s is not None and s != float('-inf'))

        actual_last_step = 0
        if all_scores:
            actual_last_step = max(all_scores.keys())

        # Now, get the best score from that actual last step
        best_last_step_score = -float('inf')
        if actual_last_step > 0:
            valid_scores_last_step = [s for s in all_scores.get(actual_last_step, []) if s is not None and s != float('-inf')]
            if valid_scores_last_step:
                best_last_step_score = max(valid_scores_last_step)
                
        best_last_step_ddg = -best_last_step_score
        # The best overall score and step are now returned directly
        best_overall_ddg = -best_overall_score
        
        # final_score = scoring_function(generated_protein)
        # final_ddg = -final_score
        
        # final_log_entry = (
        #     f"\n\n" + "="*35 + "\n    FINAL OPTIMIZED RESULT \n" + "="*35 + "\n"
        #     f"Optimized Sequence:\n{generated_protein.sequence}\n\n"
        #     f"Final Score (=-ΔΔG): {final_score:.4f}\n"
        #     f"Final FoldX ΔΔG vs. Wild-Type: {final_ddg:.4f} kcal/mol\n" + "="*35 + "\n"
        # )

        final_log_entry = (
            f"\n\n{'='*35}\n   FINAL RESULTS SUMMARY\n{'='*35}\n"
            f"\n--- Best from Final Step (Step {actual_last_step}) ---\n"
            f"Score (=-ΔΔG): {best_last_step_score:.4f}\n"
            f"FoldX ΔΔG vs. Wild-Type: {best_last_step_ddg:.4f} kcal/mol\n"
            f"\n--- Best Overall Result ---\n"
            f"Found in Step: {best_overall_step}\n"
            f"Optimized Sequence:\n{generated_protein.sequence}\n\n"
            f"Final Score (=-ΔΔG): {best_overall_score:.4f}\n"
            f"Final FoldX ΔΔG vs. Wild-Type: {best_overall_ddg:.4f} kcal/mol\n"
            f"{'='*35}\n"
        )
        
        with open(log_filepath, 'a') as f: f.write(final_log_entry)
        print(final_log_entry)
        results_dir = "../../results"
        os.makedirs(results_dir, exist_ok=True)
        plot_filename = f"ddg_history_{os.path.splitext(PDB_FILENAME)[0]}_{timestamp}.png"
        # plot_filename = f"ddg_history_{timestamp}.png"
        plot_filepath = os.path.join(results_dir, plot_filename)
        plot_ddg_history(all_scores, save_path=plot_filepath)


if __name__ == "__main__":
    main()

    