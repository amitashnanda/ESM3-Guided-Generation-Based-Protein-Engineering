from abc import ABC, abstractmethod
import os
import attr
import torch
import time
from tqdm import tqdm
from typing import Optional, Tuple, Dict
# NEW: Import Pool and cpu_count for parallel processing
from multiprocessing import Pool, cpu_count

from esm.models.esm3 import ESM3
from esm.sdk.api import (
    ESM3InferenceClient,
    ESMProtein,
    ESMProteinError,
    ESMProteinTensor,
    SamplingConfig,
    SamplingTrackConfig,
)
from esm.sdk.forge import ESM3ForgeInferenceClient
from esm.tokenization import get_esm3_model_tokenizers
from functools import partial

class GuidedDecodingScoringFunction(ABC):
    @abstractmethod
    def __call__(self, protein: ESMProtein) -> float:
        pass


class ESM3GuidedDecoding:
    """This class can be used to perform derivative-free guided decoding..."""

    def __init__(
        self,
        client: ESM3InferenceClient,
        scoring_function: GuidedDecodingScoringFunction,
    ):
        if isinstance(client, ESM3):
            self.tokenizers = client.tokenizers
        elif isinstance(client, ESM3ForgeInferenceClient):
            self.tokenizers = get_esm3_model_tokenizers(client.model)
        else:
            raise ValueError(
                "client must be an instance of ESM3 or ESM3ForgeInferenceClient"
            )

        self.client = client
        self.scoring_function = scoring_function

    def guided_generate(
        self,
        protein: ESMProtein,
        num_decoding_steps: int,
        num_samples_per_step: int,
        denoised_prediction_temperature: float = 0.0,
        track: str = "sequence",
        verbose: bool = True,
        num_workers: int = 1,
        log_file_path: Optional[str] = None
    ) -> Tuple[ESMProtein, Dict]:
        print(f"\n[START] Initial input ESMProtein:\n  sequence = {protein.sequence[:100]}...\n")
        protein_tensor = self.client.encode(protein)
        assert not isinstance(protein_tensor, ESMProteinError)

        if verbose:
            pbar = tqdm(range(num_decoding_steps), desc="Guided Generation")
        else:
            pbar = range(num_decoding_steps)

        best_overall_score = -float('inf')
        best_overall_tensor = protein_tensor
        best_overall_sequence = protein.sequence 
        best_overall_step = 0                   

        all_scores_history = {}

        for step in pbar:
            step_total_start_time = time.time()
            current_masked_positions = self.get_number_of_masked_positions(protein_tensor, track=track)
            if current_masked_positions == 0:
                print("\n[INFO] No masked positions remaining. Finishing generation.")
                break
            
            # --- MODIFIED LOGIC: PROPORTIONAL UNMASKING SCHEDULE ---
            remaining_steps = num_decoding_steps - step
            num_to_unmask = max(1, current_masked_positions // remaining_steps) if remaining_steps > 0 else current_masked_positions
            
            
            print(f"\n[STEP {step + 1}/{num_decoding_steps}] Unmasking {num_to_unmask} of {current_masked_positions} positions...")

            candidate_gen_start_time = time.time()
            candidate_tensors = [self.randomly_unmask_positions(protein_tensor, num_to_unmask, track=track) for _ in range(num_samples_per_step)]
            
            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write(f"\n--- Step {step + 1} Partial Sequences ---\n")
                    for i, tensor in enumerate(candidate_tensors):
                        partial_protein = self.client.decode(tensor)
                        f.write(f"[Candidate {i+1}] Partial Sequence: {partial_protein.sequence}\n")

            denoised_proteins = [self.predict_denoised(tensor, temperature=denoised_prediction_temperature) for tensor in candidate_tensors]
                   
            candidate_gen_duration = time.time() - candidate_gen_start_time
            print(f"[INFO] Candidate generation (GPU) finished in {candidate_gen_duration:.2f} seconds.")
        

            print(f"[INFO] Generated {len(denoised_proteins)} candidates for scoring.")

            print(f"[INFO] Scoring candidates in parallel using {num_workers} workers...")
            # step_start_time = time.time(). #previous timer
            scoring_start_time = time.time()

             # MODIFIED: Use partial to pass step and log_file_path to the scorer
            scorer_with_context = partial(
                self.scoring_function, 
                step=step + 1, 
                log_file_path=log_file_path
            )

            with Pool(processes=num_workers) as pool:
                # scores = pool.map(self.scoring_function, denoised_proteins)
                scores = pool.map(scorer_with_context, denoised_proteins)
            # step_duration = time.time() - step_start_time
            scoring_duration = time.time() - scoring_start_time
            # print(f"[INFO] Step scoring finished in {step_duration:.2f} seconds.")
            print(f"[INFO] Step scoring finished in {scoring_duration:.2f} seconds.")


            all_scores_history[step + 1] = scores

            best_score_in_step = -float('inf')
            best_tensor_in_step = None
            best_sequence_in_step = ""  #added

            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write(f"\n--- Step {step + 1} Denoised Results ---\n")
                    for i, (p, score) in enumerate(zip(denoised_proteins, scores)):
                        ddg = -score
                        ptm_str = f"pTM: {p.ptm.item():.4f}" if hasattr(p, 'ptm') and p.ptm is not None else "pTM: N/A"
                        f.write(f"\n[Candidate {i+1}]\n")
                        f.write(f"  Denoised Sequence: {p.sequence}\n")
                        f.write(f"  {ptm_str}\n")
                        f.write(f"  Score: {score:.4f} | ddG: {ddg:.3f} kcal/mol\n")
                        if score > best_score_in_step:
                            best_score_in_step = score
                            best_tensor_in_step = candidate_tensors[i]
                            best_sequence_in_step = p.sequence #added
                            
            else:
                for i, score in enumerate(scores):
                    if score > best_score_in_step:
                        best_score_in_step = score
                        best_tensor_in_step = candidate_tensors[i]
                        
            
            if best_score_in_step > best_overall_score:
                best_overall_score = best_score_in_step
                best_overall_tensor = best_tensor_in_step
                best_overall_sequence = best_sequence_in_step 
                best_overall_step = step + 1                 
            
            protein_tensor = best_tensor_in_step
            step_total_duration = time.time() - step_total_start_time

            with open(log_file_path, 'a') as f:
                f.write("\n" + "-"*25 + f" Step {step + 1} Summary " + "-"*25 + "\n")
                f.write(f"Best Candidate in Step:\n{best_sequence_in_step}\n")
                f.write(f"Best Score in Step (=-ΔΔG): {best_score_in_step:.4f}\n")
                f.write(f"This sequence will be used as the template for Step {step + 2}.\n")
                f.write("-" * 65 + "\n")

            print(f"\n[STEP {step + 1}] Best Candidate Score in Step: {best_score_in_step:.4f}")

            if verbose:
                pbar.set_description(f"Best score so far: {best_overall_score:.4f}")

            with open(log_file_path, 'a') as f:
                f.write(f"Timing for Step {step + 1}:\n")
                f.write(f"  - Candidate Generation (GPU): {candidate_gen_duration:.2f} s\n")
                f.write(f"  - Scoring (CPU/FoldX):        {scoring_duration:.2f} s\n")
                f.write(f"  - Total Step Time:              {step_total_duration:.2f} s\n")
        

            # NEW TIMER: End timer for the whole step
            # step_total_duration = time.time() - step_total_start_time
            print(f"\n[INFO] Total time for Step {step + 1}: {step_total_duration:.2f} seconds.")

        print("\n[FINAL] Performing full denoising of the best candidate to ensure completion...")
        final_complete_protein = self.predict_denoised(best_overall_tensor, temperature=0.0)
        assert not isinstance(final_complete_protein, ESMProteinError)
        
        # return final_complete_protein,all_scores_history

        return final_complete_protein, all_scores_history, best_overall_score, best_overall_step




    def reward_function(
        self,
        protein_tensor: ESMProteinTensor,
        denoised_prediction_temperature: float = 0.0,
    ) -> float:
        denoised_protein = self.predict_denoised(
            protein_tensor, temperature=denoised_prediction_temperature
        )
        print(f"[Denoised] Sequence: {denoised_protein.sequence[:60]}...")
        if hasattr(denoised_protein, "ptm"):
            print(f"[Denoised] pTM: {float(denoised_protein.ptm):.4f}")
        return self.scoring_function(denoised_protein)

    def get_number_of_masked_positions(
        self, protein_tensor: ESMProteinTensor, track: str = "sequence"
    ) -> int:
        assert isinstance(protein_tensor, ESMProteinTensor)
        track_tensor = getattr(protein_tensor, track)
        track_tokenizer = getattr(self.tokenizers, track)
        is_mask = track_tensor == track_tokenizer.mask_token_id
        return is_mask.sum().item()  # type: ignore

    def randomly_unmask_positions(
        self,
        protein_tensor: ESMProteinTensor,
        num_positions_to_unmask: int,
        temperature: float = 1.0,
        track: str = "sequence",
    ) -> ESMProteinTensor:
        track_tensor = getattr(protein_tensor, track)
        assert track_tensor is not None
        protein_tensor = attr.evolve(protein_tensor)
        setattr(protein_tensor, track, track_tensor.clone())

        track_tensor = getattr(protein_tensor, track)
        track_tokenizer = getattr(self.tokenizers, track)

        is_mask = track_tensor == track_tokenizer.mask_token_id
        num_masked_positions = is_mask.sum().item()

        if num_positions_to_unmask > num_masked_positions:
            num_positions_to_unmask = num_masked_positions  # type: ignore

        mask_indices = is_mask.nonzero(as_tuple=False)
        mask_indices = mask_indices[torch.randperm(mask_indices.size(0))]
        mask_indices = mask_indices[:num_positions_to_unmask]

        sampling_config = SamplingConfig()
        setattr(sampling_config, track, SamplingTrackConfig(temperature=temperature))

        denoised_protein_tensor_output = self.client.forward_and_sample(
            protein_tensor, sampling_configuration=sampling_config
        )
        assert not isinstance(denoised_protein_tensor_output, ESMProteinError)
        denoised_protein_tensor = denoised_protein_tensor_output.protein_tensor
        output_track_tensor = getattr(denoised_protein_tensor, track)
        assert output_track_tensor is not None
        track_tensor[mask_indices] = output_track_tensor[mask_indices]
        setattr(protein_tensor, track, track_tensor)

        return protein_tensor

    def predict_denoised(
        self, protein_tensor: ESMProteinTensor, temperature: float = 0.0
    ) -> ESMProtein:
        denoised_protein_tensor_output = self.client.forward_and_sample(
            protein_tensor,
            sampling_configuration=SamplingConfig(
                sequence=SamplingTrackConfig(temperature=temperature),
                structure=SamplingTrackConfig(temperature=temperature),
            ),
        )
        assert not isinstance(denoised_protein_tensor_output, ESMProteinError)
        denoised_protein_tensor = denoised_protein_tensor_output.protein_tensor
        denoised_protein = self.client.decode(denoised_protein_tensor)
        assert not isinstance(denoised_protein, ESMProteinError)
        return denoised_protein

    def maybe_add_default_structure_tokens(
        self, protein_tensor: ESMProteinTensor
    ) -> ESMProteinTensor:
        empty_protein_tensor = ESMProteinTensor.empty(
            len(protein_tensor) - 2,
            tokenizers=self.tokenizers,
            device=protein_tensor.device,
        )
        if protein_tensor.structure is None:
            setattr(protein_tensor, "structure", empty_protein_tensor.structure)
        else:
            print("Warning: structure already exists in protein_tensor")
        return protein_tensor
