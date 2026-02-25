from abc import ABC, abstractmethod
import gc
import os
import attr
import torch
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from typing import Optional, Tuple, Dict

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

class GuidedDecodingScoringFunction(ABC):
    @abstractmethod
    def __call__(self, protein: ESMProtein) -> float:
        pass


class ESM3GuidedDecoding:
    """This class can be used to perform derivative-free guided decoding with ESM3."""

    def __init__(
        self,
        client: ESM3InferenceClient,
        scoring_function: GuidedDecodingScoringFunction,
    ):
        # Handle DataParallel wrapped models
        if isinstance(client, torch.nn.DataParallel):
            # Access the underlying module
            actual_model = client.module
            if isinstance(actual_model, ESM3):
                self.tokenizers = actual_model.tokenizers
                self._is_data_parallel = True
            else:
                raise ValueError(
                    "DataParallel wrapped model must contain an ESM3 instance"
                )
        elif isinstance(client, ESM3):
            self.tokenizers = client.tokenizers
            self._is_data_parallel = False
        elif isinstance(client, ESM3ForgeInferenceClient):
            # All ESM3 models (small/medium/large) share the same tokenizers.
            # get_esm3_model_tokenizers only recognizes the open-small name,
            # so we call it with the default (no argument) for any Forge model.
            self.tokenizers = get_esm3_model_tokenizers()
            self._is_data_parallel = False
            self._is_forge = True
        else:
            raise ValueError(
                "client must be an instance of ESM3, DataParallel(ESM3), or ESM3ForgeInferenceClient"
            )

        self.client = client
        self.scoring_function = scoring_function
        # Set _is_forge if not already set (only True for ESM3ForgeInferenceClient)
        if not hasattr(self, '_is_forge'):
            self._is_forge = False

    @property
    def _model(self):
        """Get the actual model, unwrapping DataParallel if needed."""
        if self._is_data_parallel:
            return self.client.module
        return self.client

    def guided_generate(
        self,
        protein: ESMProtein,
        num_decoding_steps: int,
        num_samples_per_step: int,
        denoised_prediction_temperature: float = 0.0,
        track: str = "sequence",
        verbose: bool = True,
        log_file_path: Optional[str] = None,
        # num_boltz_workers: int = 1
    ) -> Tuple[ESMProtein, Dict]:
        print(f"\n[START] Initial input ESMProtein:\n  sequence = {protein.sequence[:100]}...\n")
        protein_tensor = self._model.encode(protein)
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
            
            # --- PROPORTIONAL UNMASKING SCHEDULE ---
            remaining_steps = num_decoding_steps - step
            num_to_unmask = max(1, current_masked_positions // remaining_steps) if remaining_steps > 0 else current_masked_positions
            
            
            print(f"\n[STEP {step + 1}/{num_decoding_steps}] Unmasking {num_to_unmask} of {current_masked_positions} positions...")

            candidate_gen_start_time = time.time()
            
            # ================================================================
            # FORGE API: Parallelize candidate generation using ThreadPoolExecutor
            # Each forward_and_sample call is an independent HTTP request, so
            # running them concurrently gives ~Nx speedup with N threads.
            # LOCAL MODEL: Keep sequential since DataParallel already uses all
            # GPUs within each forward pass.
            # ================================================================
            if self._is_forge:
                # --- PARALLEL candidate generation (Forge API) ---
                print(f"[INFO] Generating {num_samples_per_step} candidates in parallel (Forge API)...")
                with ThreadPoolExecutor(max_workers=num_samples_per_step) as executor:
                    unmask_futures = [executor.submit(self.randomly_unmask_positions, protein_tensor, num_to_unmask, 1.0, track) 
                                      for _ in range(num_samples_per_step)]
                    candidate_tensors = [f.result() for f in unmask_futures]
            else:
                # --- SEQUENTIAL candidate generation (Local model) ---
                # # OLD SEQUENTIAL CODE (commented out for reference):
                # candidate_tensors = [self.randomly_unmask_positions(protein_tensor, num_to_unmask, track=track) for _ in range(num_samples_per_step)]
                candidate_tensors = [self.randomly_unmask_positions(protein_tensor, num_to_unmask, track=track) for _ in range(num_samples_per_step)]
            
            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write(f"\n--- Step {step + 1} Partial Sequences ---\n")
                    for i, tensor in enumerate(candidate_tensors):
                        partial_protein = self._model.decode(tensor)
                        f.write(f"[Candidate {i+1}] Partial Sequence: {partial_protein.sequence}\n")

            # Create the denoised proteins
            if self._is_forge:
                # --- PARALLEL denoised prediction (Forge API) ---
                print(f"[INFO] Denoising {len(candidate_tensors)} candidates in parallel (Forge API)...")
                with ThreadPoolExecutor(max_workers=len(candidate_tensors)) as executor:
                    denoise_futures = [executor.submit(self.predict_denoised, tensor, denoised_prediction_temperature) 
                                       for tensor in candidate_tensors]
                    raw_denoised_proteins = [f.result() for f in denoise_futures]
            else:
                # --- SEQUENTIAL denoised prediction (Local model) ---
                # # OLD SEQUENTIAL CODE (commented out for reference):
                # raw_denoised_proteins = [self.predict_denoised(tensor, temperature=denoised_prediction_temperature) for tensor in candidate_tensors]
                raw_denoised_proteins = [self.predict_denoised(tensor, temperature=denoised_prediction_temperature) for tensor in candidate_tensors]

            denoised_proteins = []
            for protein_item in raw_denoised_proteins:
                detached_protein = ESMProtein(sequence=protein_item.sequence)
                for attr_name, attr_value in vars(protein_item).items():
                    if isinstance(attr_value, torch.Tensor):
                        setattr(detached_protein, attr_name, attr_value.detach().clone()) 
                    elif attr_name != 'sequence': 
                         setattr(detached_protein, attr_name, attr_value) 
                denoised_proteins.append(detached_protein)
                   
            candidate_gen_duration = time.time() - candidate_gen_start_time
            print(f"[INFO] Candidate generation (GPU) finished in {candidate_gen_duration:.2f} seconds.")
        
            # Clear GPU memory before Boltz scoring to prevent contention
            gc.collect()
            torch.cuda.empty_cache()
            
            print(f"[INFO] Generated {len(denoised_proteins)} candidates for scoring.")
            
            scoring_start_time = time.time()

            # ================================================================
            # BATCHED SCORING: Score all candidates in ONE Boltz subprocess call
            # This loads the Boltz model once per step instead of once per candidate
            # (~3x faster than per-candidate scoring)
            # ================================================================
            print(f"[INFO] Scoring all {len(denoised_proteins)} candidates in one batched Boltz call...")
            
            # Use score_batch if available (BoltzScorer), otherwise fall back to individual
            if hasattr(self.scoring_function, 'score_batch'):
                batch_results = self.scoring_function.score_batch(
                    proteins=denoised_proteins,
                    step=step + 1,
                    log_file_path=log_file_path,
                    # num_devices auto-detected: uses all GPUs on node
                )
                scores = [r[0] for r in batch_results]
                affinity_details = [r[1] for r in batch_results]
                for idx, (score, _) in enumerate(batch_results):
                    if score == float('-inf'):
                        print(f"  Scored candidate {idx+1}/{len(denoised_proteins)}: -inf")
                    else:
                        print(f"  Scored candidate {idx+1}/{len(denoised_proteins)}: {score:.4f}")
            else:
                # Fallback: score individually if scoring function doesn't support batching
                print(f"[INFO] Scoring function does not support batching, falling back to sequential...")
                scores = []
                affinity_details = []
                for idx, protein_item in enumerate(denoised_proteins):
                    print(f"  Scoring candidate {idx+1}/{len(denoised_proteins)}...")
                    score, details = self.scoring_function(protein_item, step=step + 1, log_file_path=log_file_path)
                    scores.append(score)
                    affinity_details.append(details)

            # ================================================================
            # OLD PER-CANDIDATE SCORING (COMMENTED OUT)
            # This spawned a separate subprocess per candidate, reloading the
            # Boltz model each time. Replaced by batched scoring above.
            # ================================================================
            # # Score candidates - parallel or sequential based on num_boltz_workers
            # # Detect available GPUs for round-robin assignment
            # num_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 1
            # 
            # if num_boltz_workers > 1:
            #     print(f"[INFO] Scoring candidates in parallel ({num_boltz_workers} workers across {num_gpus} GPUs)...")
            #     scores = [None] * len(denoised_proteins)
            #     affinity_details = [None] * len(denoised_proteins)
            #     
            #     def score_single(idx_protein_gpu):
            #         idx, protein_item, assigned_gpu = idx_protein_gpu
            #         return idx, self.scoring_function(protein_item, step=step + 1, log_file_path=log_file_path, gpu_id=assigned_gpu)
            #     
            #     with ThreadPoolExecutor(max_workers=num_boltz_workers) as executor:
            #         # Round-robin GPU assignment: candidate 0 → GPU 0, candidate 1 → GPU 1, etc.
            #         futures = {executor.submit(score_single, (i, p, i % num_gpus)): i 
            #                    for i, p in enumerate(denoised_proteins)}
            #         for future in as_completed(futures):
            #             idx, result = future.result()
            #             score, details = result
            #             scores[idx] = score
            #             affinity_details[idx] = details
            #             print(f"  Scored candidate {idx+1}/{len(denoised_proteins)} (GPU {idx % num_gpus}): {score:.4f}")
            # else:
            #     # Sequential scoring (original behavior)
            #     print(f"[INFO] Scoring candidates sequentially (Boltz uses GPU)...")
            #     scores = []
            #     affinity_details = []
            #     for idx, protein_item in enumerate(denoised_proteins):
            #         print(f"  Scoring candidate {idx+1}/{len(denoised_proteins)}...")
            #         score, details = self.scoring_function(protein_item, step=step + 1, log_file_path=log_file_path)
            #         scores.append(score)
            #         affinity_details.append(details)
            # ================================================================

            scoring_duration = time.time() - scoring_start_time
            print(f"[INFO] Step scoring finished in {scoring_duration:.2f} seconds.")

            all_scores_history[step + 1] = scores

            best_score_in_step = -float('inf')
            best_tensor_in_step = None
            best_sequence_in_step = ""

            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write(f"\n--- Step {step + 1} Denoised Results ---\n")
                    for i, (p, score) in enumerate(zip(denoised_proteins, scores)):
                        ptm_str = f"pTM: {p.ptm.item():.4f}" if hasattr(p, 'ptm') and p.ptm is not None else "pTM: N/A"
                        details = affinity_details[i] if i < len(affinity_details) else None
                        
                        f.write(f"\n[Candidate {i+1}]\n")
                        f.write(f"  Denoised Sequence: {p.sequence}\n")
                        f.write(f"  {ptm_str}\n")
                        
                        # Print raw Boltz value and negated score side-by-side
                        if details is not None:
                            raw_val = details.get('affinity_pred_value', 'N/A')
                            f.write(f"  Score: {score:.4f} | Raw Boltz affinity_pred_value (log10 IC50): {raw_val}\n")
                            f.write(f"  --- Boltz Affinity Details ---\n")
                            f.write(f"    affinity_pred_value (ensemble):  {details.get('affinity_pred_value', 'N/A')}\n")
                            f.write(f"    affinity_probability_binary:     {details.get('affinity_probability_binary', 'N/A')}\n")
                            f.write(f"    affinity_pred_value1 (model 1):  {details.get('affinity_pred_value1', 'N/A')}\n")
                            f.write(f"    affinity_probability_binary1:    {details.get('affinity_probability_binary1', 'N/A')}\n")
                            f.write(f"    affinity_pred_value2 (model 2):  {details.get('affinity_pred_value2', 'N/A')}\n")
                            f.write(f"    affinity_probability_binary2:    {details.get('affinity_probability_binary2', 'N/A')}\n")
                        else:
                            f.write(f"  Score: {score:.4f} | Boltz prediction failed\n")
                        
                        if score > best_score_in_step:
                            best_score_in_step = score
                            best_tensor_in_step = candidate_tensors[i]
                            best_sequence_in_step = p.sequence
                            
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
            
            # If all candidates scored -inf (e.g., Boltz failed), keep the
            # previous step's tensor so we can continue to the next step
            if best_tensor_in_step is None:
                print(f"[WARNING] All candidates in Step {step + 1} scored -inf. "
                      f"Retaining previous step's best protein for next step.")
                # protein_tensor stays unchanged from the previous iteration
            else:
                protein_tensor = best_tensor_in_step
            step_total_duration = time.time() - step_total_start_time

            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write("\n" + "-"*25 + f" Step {step + 1} Summary " + "-"*25 + "\n")
                    f.write(f"Best Candidate in Step:\n{best_sequence_in_step}\n")
                    f.write(f"Best Score in Step (=-affinity_pred_value): {best_score_in_step:.4f}\n")
                    f.write(f"Raw Boltz affinity_pred_value (log10 IC50): {-best_score_in_step:.4f}\n")
                    f.write(f"This sequence will be used as the template for Step {step + 2}.\n")
                    f.write("-" * 65 + "\n")

            print(f"\n[STEP {step + 1}] Best Candidate Score in Step: {best_score_in_step:.4f}")

            if verbose:
                pbar.set_description(f"Best score so far: {best_overall_score:.4f}")

            if log_file_path:
                with open(log_file_path, 'a') as f:
                    f.write(f"Timing for Step {step + 1}:\n")
                    f.write(f"  - Candidate Generation (GPU): {candidate_gen_duration:.2f} s\n")
                    f.write(f"  - Scoring (GPU/Boltz):        {scoring_duration:.2f} s\n")
                    f.write(f"  - Total Step Time:            {step_total_duration:.2f} s\n")
        
            print(f"\n[INFO] Total time for Step {step + 1}: {step_total_duration:.2f} seconds.")

        print("\n[FINAL] Performing full denoising of the best candidate to ensure completion...")
        final_complete_protein = self.predict_denoised(best_overall_tensor, temperature=0.0)
        assert not isinstance(final_complete_protein, ESMProteinError)
        
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

        denoised_protein_tensor_output = self._model.forward_and_sample(
            protein_tensor, sampling_configuration=sampling_config
        )
        assert not isinstance(denoised_protein_tensor_output, ESMProteinError)
        denoised_protein_tensor = denoised_protein_tensor_output.protein_tensor
        output_track_tensor = getattr(denoised_protein_tensor, track).long()
        assert output_track_tensor is not None
        track_tensor[mask_indices] = output_track_tensor[mask_indices].long()
        setattr(protein_tensor, track, track_tensor)

        return protein_tensor

    def predict_denoised(
        self, protein_tensor: ESMProteinTensor, temperature: float = 0.0
    ) -> ESMProtein:
        denoised_protein_tensor_output = self._model.forward_and_sample(
            protein_tensor,
            sampling_configuration=SamplingConfig(
                sequence=SamplingTrackConfig(temperature=temperature),
                structure=SamplingTrackConfig(temperature=temperature),
            ),
        )
        assert not isinstance(denoised_protein_tensor_output, ESMProteinError)
        denoised_protein_tensor = denoised_protein_tensor_output.protein_tensor
        denoised_protein = self._model.decode(denoised_protein_tensor)
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



