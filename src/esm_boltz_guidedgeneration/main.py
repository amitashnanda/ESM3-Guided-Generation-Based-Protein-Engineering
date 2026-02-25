"""
ESM3-Boltz Guided Generation - Main Script

This script performs ESM3 guided generation to design protein sequences
with enhanced binding affinity to a specified ligand using Boltz-2 scoring.

Usage:
    python main.py --smiles "CCO" --seq_length 256 --masking_percentage 0.4 \\
                   --num_decoding_steps 32 --num_samples_per_step 20

Author: Amitash Nanda
Last modified by Anna Su
Date: 2026-02-20
"""

import os
import random
import time
import argparse
import warnings

import torch
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein
from esm.utils.structure.protein_chain import ProteinChain

# Import Forge API client (optional - only needed for remote models)
try:
    from esm.sdk.forge import ESM3ForgeInferenceClient
    FORGE_AVAILABLE = True
except ImportError:
    FORGE_AVAILABLE = False

# Import custom modules (must be in same directory or PYTHONPATH)
import sys
sys.path.insert(0, os.path.dirname(__file__))
from guided_generation import ESM3GuidedDecoding
from boltz_scoring_utils import BoltzScorer, plot_affinity_history


warnings.simplefilter(action='ignore', category=FutureWarning)


# --- Default Settings ---
TIMEOUT_SEC = 7200
CLEANUP_TMP = True
VERBOSE_BOLTZ = False

DEFAULT_LOG_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "logs", "boltz")
)

DEFAULT_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "results", "boltz")
)


def get_masked_sequence(wildtype_pdb: str, wildtype_chain: str, masking_percentage: float) -> str:
    """
    Load a wildtype sequence from PDB and create a masked version.
    
    Args:
        wildtype_pdb: PDB ID (e.g., "1ABC")
        wildtype_chain: Chain ID (e.g., "A")
        masking_percentage: Fraction of residues to mask (0.0 to 1.0)
    
    Returns:
        Masked sequence string
    """
    print("[INFO] Loading wildtype from PDB: {}, Chain: {}".format(wildtype_pdb, wildtype_chain))
    
    wildtype_protein = ESMProtein.from_protein_chain(
        ProteinChain.from_rcsb(wildtype_pdb, chain_id=wildtype_chain)
    )
    
    wt_sequence = wildtype_protein.sequence
    maskable_indices = list(range(len(wt_sequence)))
    num_to_mask = int(len(maskable_indices) * masking_percentage)
    indices_to_mask = random.sample(maskable_indices, num_to_mask)
    
    refinement_template_list = list(wt_sequence)
    for i in indices_to_mask:
        refinement_template_list[i] = '_'
    
    masked_seq = "".join(refinement_template_list)
    
    print("[INFO] Wildtype sequence length: {}".format(len(wt_sequence)))
    print("[INFO] Masked {} positions ({:.1f}%)".format(num_to_mask, masking_percentage * 100))
    print("[INFO] Masked sequence: {}...".format(masked_seq[:80]))
    
    return masked_seq


def main():
    parser = argparse.ArgumentParser(
        description="Run ESM3 Guided Generation with Boltz-1 binding affinity scoring.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "--smiles",
        type=str,
        required=True,
        help="SMILES string of the target ligand for binding affinity prediction."
    )
    
    # ESM3 Model selection
    parser.add_argument(
        "--model",
        type=str,
        choices=["local", "esm3-medium-2024-08", "esm3-large"],
        default="local",
        help="ESM3 model to use. 'local' uses esm3-open-small on GPU, others use Forge API. Default: local"
    )
    parser.add_argument(
        "--forge_token",
        type=str,
        default=None,
        help="Forge API token (required for esm3-medium/large). Can also set ESM_FORGE_TOKEN env var."
    )
    
    # Sequence initialization options
    sequence_group = parser.add_mutually_exclusive_group(required=True)
    sequence_group.add_argument(
        "--seq_length",
        type=int,
        help="Length of de novo protein sequence to generate (fully masked start)."
    )
    sequence_group.add_argument(
        "--wildtype",
        nargs=2,
        metavar=("PDB_ID", "CHAIN_ID"),
        help="Start from wildtype PDB structure (e.g., --wildtype 1ABC A)."
    )
    
    # Generation parameters
    parser.add_argument(
        "--masking_percentage",
        type=float,
        default=0.40,
        help="Percentage of residues to mask (only used with --wildtype). Default: 0.40"
    )
    parser.add_argument(
        "--num_decoding_steps",
        type=int,
        default=32,
        help="Number of generation steps. Default: 32"
    )
    parser.add_argument(
        "--num_samples_per_step",
        type=int,
        default=20,
        help="Number of candidates to generate per step. Default: 20"
    )
    
    # Output options
    parser.add_argument(
        "--log_file",
        type=str,
        default=None,
        help="Path to generation log file. Defaults to timestamped file in logs/boltz/"
    )
    parser.add_argument(
        "--run_dir",
        type=str,
        default=None,
        help="Run directory name (e.g., run_20260220_204005). When set, logs and results are placed under this directory to keep .out/.err/.txt together."
    )
    parser.add_argument(
        "--cache_dir",
        type=str,
        default=None,
        help="Directory for caching Boltz results. Default: ./boltz_cache"
    )
    parser.add_argument(
        "--score_type",
        type=str,
        choices=["affinity", "binary"],
        default="affinity",
        help="Score type: 'affinity' for lead optimization (log10 IC50), 'binary' for hit discovery (probability). Default: affinity"
    )
    parser.add_argument(
        "--use_msa_server",
        action="store_true",
        help="Use Boltz MSA server (slower but more accurate). Computes MSA per candidate sequence."
    )
    # parser.add_argument(
    #     "--num_boltz_workers",
    #     type=int,
    #     default=1,
    #     help="Number of parallel Boltz workers for scoring. Default: 1 (sequential). Set higher for parallel scoring."
    # )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.wildtype and not (0.0 < args.masking_percentage < 1.0):
        parser.error("--masking_percentage must be between 0.0 and 1.0")
    
    if args.seq_length and args.seq_length < 1:
        parser.error("--seq_length must be positive")
    
    # Setup
    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    
    # Create run directory for this execution
    # If --run_dir is provided (e.g., from guided_generation.sh), use it so that
    # .txt log files end up in the same directory as .out/.err files.
    # Fallback is set after wildtype parsing to include PDB ID in the name.
    run_dir_name = args.run_dir  # May be None; finalized below
    
    # Setup cache directory
    if args.cache_dir:
        cache_dir = args.cache_dir
    else:
        cache_dir = os.path.join(os.getcwd(), "boltz_cache")
    os.makedirs(cache_dir, exist_ok=True)
    
    # --- 1. Create Starting Sequence ---
    if args.wildtype:
        wildtype_pdb, wildtype_chain = args.wildtype
        masked_seq = get_masked_sequence(wildtype_pdb, wildtype_chain, args.masking_percentage)
        seq_identifier = "{}_{}".format(wildtype_pdb, wildtype_chain)
    else:
        # De novo generation - fully masked sequence
        masked_seq = "_" * args.seq_length
        seq_identifier = "denovo_{}".format(args.seq_length)
        print("[INFO] Starting de novo generation with sequence length: {}".format(args.seq_length))
    
    # Finalize run directory name (include PDB ID for easier identification)
    if not run_dir_name:
        if args.wildtype:
            run_dir_name = "run_{}_{}".format(wildtype_pdb, timestamp)
        else:
            run_dir_name = "run_denovo_{}".format(timestamp)
    
    starting_protein = ESMProtein(sequence=masked_seq)
    
    # --- 2. Setup Log File ---
    if args.log_file:
        log_filepath = args.log_file
    else:
        # Create run-specific log directory
        run_log_dir = os.path.join(DEFAULT_LOG_DIR, run_dir_name)
        os.makedirs(run_log_dir, exist_ok=True)
        mask_str = "{}".format(int(args.masking_percentage * 100)) if args.wildtype else "100"
        # Use the original timestamp (not run_dir) for the filename so it's unique per job
        log_filename = "{}_{}_mask{}_steps{}_{}.txt".format(seq_identifier, args.model, mask_str, args.num_decoding_steps, timestamp)
        log_filepath = os.path.join(run_log_dir, log_filename)
    
    log_dirname = os.path.dirname(log_filepath)
    if log_dirname:
        os.makedirs(log_dirname, exist_ok=True)
    
    print("[INFO] Log file: {}".format(log_filepath))
    
    # Write header to log file
    header_text = (
        "======================================================================\n"
        "ESM3-Boltz Guided Generation Log\n"
        "Run Timestamp: {}\n".format(timestamp) +
        "======================================================================\n\n"
        "--- RUN PARAMETERS ---\n"
        "ESM3 Model: {}\n".format(args.model) +
        "Target Ligand SMILES: {}\n".format(args.smiles)
    )
    
    if args.wildtype:
        header_text += (
            "Starting from: PDB {}, Chain {}\n".format(wildtype_pdb, wildtype_chain) +
            "Masking Percentage: {:.1f}%\n".format(args.masking_percentage * 100)
        )
    else:
        header_text += (
            "De novo generation\n"
            "Sequence Length: {}\n".format(args.seq_length)
        )
    
    header_text += (
        "Decoding Steps: {}\n".format(args.num_decoding_steps) +
        "Samples per Step: {}\n".format(args.num_samples_per_step) +
        "Boltz Cache Dir: {}\n".format(cache_dir) +
        # "Boltz Workers: {}\n".format(args.num_boltz_workers) +
        "Score Type: {}\n".format(args.score_type) +
        "MSA Server: {}\n".format("Enabled (per-candidate)" if args.use_msa_server else "Disabled") +
        "Boltz will use GPU for batched parallel scoring\n\n"
        "--- SEQUENCE SETUP ---\n"
        "Starting Masked Sequence:\n{}\n".format(masked_seq) +
        "======================================================================\n"
    )
    
    with open(log_filepath, 'w') as f:
        f.write(header_text)
    print(header_text)
    
    # --- 3. Initialize Model ---
    if args.model == "local":
        # Use local ESM3-open-small model
        device = "cuda" if torch.cuda.is_available() else "cpu"
        print("\nLoading local ESM3-open-small model to {}...".format(device))
        
        num_gpus = torch.cuda.device_count()
        print("Available GPUs: {}".format(num_gpus))
        
        model = ESM3.from_pretrained().float()
        
        # Use DataParallel for multi-GPU
        if num_gpus > 1:
            print("Using DataParallel across {} GPUs: {}".format(num_gpus, [torch.cuda.get_device_name(i) for i in range(num_gpus)]))
            model = torch.nn.DataParallel(model)
            model = model.to(device)
        else:
            model = model.to(device)
            if num_gpus == 1:
                print("Using single GPU: {}".format(torch.cuda.get_device_name(0)))
            else:
                print("Using CPU")
    else:
        # Use Forge API for esm3-medium or esm3-large
        if not FORGE_AVAILABLE:
            raise RuntimeError(
                "Forge API client not available. Install with: pip install esm[forge]\n"
                "Or use --model local to run locally."
            )
        
        # Get token from argument or environment variable
        forge_token = args.forge_token or os.environ.get("ESM_FORGE_TOKEN")
        if forge_token:
            forge_token = forge_token.strip()  # Remove whitespace to prevent "Illegal header value" errors
        if not forge_token:
            raise ValueError(
                "Forge API token required for {}. Provide via --forge_token or set ESM_FORGE_TOKEN env var.\n"
                "Get your token at: https://forge.evolutionaryscale.ai".format(args.model)
            )
        
        print("\nConnecting to Forge API for {}...".format(args.model))
        model = ESM3ForgeInferenceClient(
            model=args.model,
            url="https://forge.evolutionaryscale.ai",
            token=forge_token
        )
        print("Connected to Forge API successfully!")
    
    # --- 4. Initialize Boltz Scorer ---
    print("\nInitializing Boltz scorer...")
    print("Cache directory: {}".format(cache_dir))
    
    msa_info = "per-candidate" if args.use_msa_server else "disabled"
    print("Boltz will use GPU for batched parallel scoring (MSA: {})".format(msa_info))
    
    scoring_function = BoltzScorer(
        smiles=args.smiles,
        cache_dir=cache_dir,
        timeout_sec=TIMEOUT_SEC,
        cleanup_tmp=CLEANUP_TMP,
        verbose_boltz=VERBOSE_BOLTZ,
        score_type=args.score_type,
        use_msa_server=args.use_msa_server,
    )
    
    # --- 5. Run Guided Generation ---
    guided_decoding = ESM3GuidedDecoding(client=model, scoring_function=scoring_function)
    
    print("\n======================================================================")
    print("Starting ESM3-Boltz Guided Generation")
    print("======================================================================\n")
    
    generated_protein, all_scores, best_overall_score, best_overall_step = guided_decoding.guided_generate(
        protein=starting_protein,
        num_decoding_steps=args.num_decoding_steps,
        num_samples_per_step=args.num_samples_per_step,
        track="sequence",
        log_file_path=log_filepath,
        # num_boltz_workers=args.num_boltz_workers
    )
    
    # --- 6. Display and Save Results ---
    if generated_protein:
        # Get best score from last step
        actual_last_step = max(all_scores.keys()) if all_scores else 0
        best_last_step_score = -float('inf')
        
        if actual_last_step > 0:
            valid_scores_last_step = [
                s for s in all_scores.get(actual_last_step, [])
                if s is not None and s != float('-inf')
            ]
            if valid_scores_last_step:
                best_last_step_score = max(valid_scores_last_step)
        
        # Write final summary
        final_log_entry = (
            "\n\n======================================================================\n"
            "   FINAL RESULTS SUMMARY\n"
            "======================================================================\n\n"
            "--- Best from Final Step (Step {}) ---\n".format(actual_last_step) +
            "Score (=-affinity_pred_value): {:.4f}\n".format(best_last_step_score) +
            "Boltz affinity_pred_value (log10 IC50): {:.4f}\n\n".format(-best_last_step_score) +
            "--- Best Overall Result ---\n"
            "Found in Step: {}\n".format(best_overall_step) +
            "Optimized Sequence:\n{}\n\n".format(generated_protein.sequence) +
            "Final Score (=-affinity_pred_value): {:.4f}\n".format(best_overall_score) +
            "Final Boltz affinity_pred_value (log10 IC50): {:.4f}\n".format(-best_overall_score) +
            "======================================================================\n"
        )
        
        with open(log_filepath, 'a') as f:
            f.write(final_log_entry)
        print(final_log_entry)
        
        # Save plots in run-specific directory
        run_results_dir = os.path.join(DEFAULT_RESULTS_DIR, run_dir_name)
        os.makedirs(run_results_dir, exist_ok=True)
        
        # Use new plotting function from notebook
        plot_affinity_history(log_filepath, run_results_dir, seq_identifier)
        
        print("\n✓ Generation complete!")
        print("✓ Log saved to: {}".format(log_filepath))
        print("✓ Plots saved to: {}".format(run_results_dir))
    else:
        print("[ERROR] Generation failed to produce a valid protein.")


if __name__ == "__main__":
    main()


# Example usage:
# De novo generation (local model):
#   python main.py --smiles "CCO" --seq_length 256 --num_decoding_steps 32 --num_samples_per_step 20
#   -> Log: logs/boltz/run_denovo_<timestamp>/denovo_256_local_mask100_steps32_<timestamp>.txt
#
# From wildtype (local model):
#   python main.py --smiles "NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O" \
#       --wildtype 1RNT A --masking_percentage 0.4 --num_decoding_steps 32 --num_samples_per_step 10
#   -> Log: logs/boltz/run_1RNT_<timestamp>/1RNT_A_local_mask40_steps32_<timestamp>.txt
#
# Using Forge API with MSA server (per-candidate MSA):
#   python main.py --model esm3-medium-2024-08 --forge_token YOUR_TOKEN \
#       --smiles "NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O" \
#       --wildtype 1RNT A --masking_percentage 0.4 --num_decoding_steps 32 --num_samples_per_step 10 \
#       --cache_dir "./boltz_cache/1rnt" --score_type affinity --use_msa_server
#   -> Token whitespace is auto-stripped to prevent "Illegal header value" errors.
#
# Binary scoring (hit discovery):
#   python main.py --model esm3-medium-2024-08 --forge_token YOUR_TOKEN \
#       --smiles "CCO" --wildtype 1DB1 A --score_type binary --num_decoding_steps 32 \
#       --num_samples_per_step 10 --cache_dir "./boltz_cache/1db1"
#
# With explicit run directory (used internally by guided_generation.sh to co-locate logs):
#   python main.py --smiles "CCO" --wildtype 1RNT A --run_dir "run_1RNT_20260220_204005"
#   -> All logs go to logs/boltz/run_1RNT_20260220_204005/ instead of a new timestamped dir.