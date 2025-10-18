"""This script uses ESM-3 guided generation to design protein sequences
with enhanced binding affinity to a specified ligand using Boltz scoring.

Warning: 
The script has not been tested since we currently do not have a single environment 
for boltz, esm3, and foldx due to dependency conflicts.

Author: Anna Su <anna.su@yale.edu>
Date: 2025-10-17
"""

import argparse as ap
import random
import argparse as ap
from datetime import datetime
from pathlib import Path

from esm.sdk.api import ESMProtein
from esm.sdk.experimental import ESM3GuidedDecoding, GuidedDecodingScoringFunction
from esm.utils.structure.protein_chain import ProteinChain
from esm.sdk import client

from boltz_scoring_utils import get_boltz_affinity

from esm.models.esm3 import ESM3
model = ESM3.from_pretrained().to("cuda")

# Binding affinity scoring function using Boltz
class AffinityScoringFunction(GuidedDecodingScoringFunction):
    def __init__(self, out_dir, smiles):
             super().__init__()
             self.out_dir = out_dir
             self.smiles = smiles

    def __call__(self, protein: ESMProtein) -> float:
        score = self.affinity(protein)
        return score

    def affinity(self, protein: ESMProtein) -> float:
        seq = protein.sequence
        try:
            affinity = get_boltz_affinity(seq, self.smiles, self.out_dir)
            print(f"******Protein seq for binding affinity {affinity}: {seq}*******")
            return affinity
        except Exception as e:
            print(f"Boltz binding affinity prediction failed: {e}")
            return 0.0

def get_masked_sequence(wildtype, masking_percentage):
    wildtype_pdb = wildtype[0]
    wildtype_chain = wildtype[1]
    wildtype = ESMProtein.from_protein_chain(
            ProteinChain.from_rcsb(wildtype_pdb, chain_id=wildtype_chain)
        )
    maskable_indices = list(range(len(wildtype.sequence)))
    num_to_mask = int(len(maskable_indices) * masking_percentage)
    indices_to_mask = random.sample(maskable_indices, num_to_mask)
    refinement_template_list = list(wildtype.sequence)

    for i in indices_to_mask: 
        refinement_template_list[i] = '_'
    masked_seq = "".join(refinement_template_list)
    print(f"Masked sequence: {masked_seq}")
    return masked_seq

def run_guided_generate(out_dir, seq_len, smiles, wildtype, num_decoding_steps, num_samples_per_step, masking_percentage):

    if wildtype is not None:
        # Start from a wildtype sequence
        masked_seq = get_masked_sequence(wildtype, masking_percentage)
        starting_protein = ESMProtein(sequence=masked_seq)
    else:
        # Start from a fully masked protein
        starting_protein = ESMProtein(sequence="_" * seq_len)

    scoring_fn = AffinityScoringFunction(out_dir, smiles)
    affinity_guided_decoding = ESM3GuidedDecoding(client=model, scoring_function=scoring_fn)

    # Call guided_generate
    generated_protein = affinity_guided_decoding.guided_generate(
        protein=starting_protein,
        num_decoding_steps=num_decoding_steps,
        num_samples_per_step=num_samples_per_step
    )

    return generated_protein

def main():
    parser = ap.ArgumentParser()
    parser.add_argument("--out_dir", type=str, help="Output directory to save generated PDBs")
    parser.add_argument("--seq_length", default=256, type=int, help="Sequence Length for the generated protein")
    parser.add_argument("--smiles", default=None, type=str, help="SMILES representation of the ligand")
    parser.add_argument("--wildtype", nargs=2, metavar=("PDB_ID", "CHAIN_ID"), default=None, help="Provide PDB ID and chain (e.g. --wildtype 1uak A) to initiate the guided generation")
    parser.add_argument("--num_decoding_steps", default=64, type=int, help="Number of steps in the guided generation")
    parser.add_argument("--num_samples_per_step", default=10, type=int, help="Number of sample proteins generated in each step")
    parser.add_argument("--masking_percentage", default=0.4, type=float, help="Percentage of the sequence to mask")

    args = parser.parse_args()

    # Build a default out_dir if none provided: runs/run_YYYYmmdd_HHMMSS
    if not args.out_dir:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_dir = Path("runs") / f"run_{ts}"
        out_dir = str(default_dir)
    else:
        out_dir = args.out_dir

    # Ensure output directory exists
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out_dir}")

    # Call guided generation
    run_guided_generate(
        out_dir,
        args.seq_length,
        args.smiles,
        args.wildtype,
        args.num_decoding_steps,
        args.num_samples_per_step,
        args.masking_percentage
    )

if __name__=='__main__':
    main()

# Example usage:
# python guided_generation.py --seq_length 400 --smiles "CCO" --wildtype 2TRX A