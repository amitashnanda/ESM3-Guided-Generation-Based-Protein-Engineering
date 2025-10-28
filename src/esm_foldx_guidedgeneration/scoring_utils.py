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
from guided_generation import ESM3GuidedDecoding, GuidedDecodingScoringFunction
from multiprocessing import Pool, cpu_count, Manager
from functools import partial
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from multiprocessing import set_start_method
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


AA3_TO_AA1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def parse_pdb_chain_sequence_with_mapping(pdb_path: str, chain_id: str) -> Tuple[str, List[str]]:
    seq, mapping, seen = [], [], set()
    with open(pdb_path, "r") as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[21].strip() != chain_id: continue
            resname = line[17:20].strip().upper()
            if resname not in AA3_TO_AA1: continue
            resseq, icode = line[22:26].strip(), line[26].strip()
            key = (resseq, icode)
            if key in seen: continue
            seen.add(key)
            seq.append(AA3_TO_AA1[resname])
            mapping.append(resseq + (icode if icode else ""))
    return "".join(seq), mapping

def foldx_repair_pdb(foldx_exec: str, workdir: str, pdb_filename: str,timeout_sec: int) -> str:
    print(f"[INFO] Running FoldX RepairPDB on {pdb_filename}...")
    cmd = [foldx_exec, "--command=RepairPDB", f"--pdb={pdb_filename}"]
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=timeout_sec)
    if proc.returncode != 0: raise RuntimeError(f"RepairPDB failed. FoldX stderr:\n{proc.stderr}")
    repaired_path = os.path.join(workdir, f"{os.path.splitext(pdb_filename)[0]}_Repair.pdb")
    if not os.path.isfile(repaired_path): raise FileNotFoundError(f"RepairPDB finished but output file not found: {repaired_path}")
    print(f"[OK] PDB repaired. Output: {repaired_path}")
    return repaired_path

def get_cache_path(seq: str, cache_dir: str) -> str:
    h = hashlib.sha256(seq.encode()).hexdigest()
    return os.path.join(cache_dir, f"ddg_{h[:16]}.joblib")

def parse_ddg_from_fxout(fxout_path: str) -> Optional[float]:
    try:
        with open(fxout_path, "r") as f: lines = [line.strip() for line in f if line.strip()]
        header_line_index, header = -1, []
        for i, line in enumerate(lines):
            if "pdb" in line.lower() and ("total energy" in line.lower() or "ddg" in line.lower()):
                header_line_index, header = i, line.split('\t')
                break
        if header_line_index == -1: return None
        data_line_index = header_line_index + 1
        if len(lines) <= data_line_index: return None
        data_line, possible_colnames = lines[data_line_index].split('\t'), ["ddg", "total energy", "energydiff", "interaction energy"]
        energy_col_idx = -1
        for i, col_name in enumerate(header):
            if col_name.lower().strip().replace("_", " ") in possible_colnames:
                energy_col_idx = i
                break
        if energy_col_idx == -1 or len(data_line) <= energy_col_idx: return None
        return float(data_line[energy_col_idx])
    except Exception as e:
        print(f"[WARN] Could not parse ddG from {fxout_path}: {e}")
        return None

def get_foldx_mutation_string(candidate_seq: str, wt_seq: str, chain_id: str, seq_to_pdb_map: List[str]) -> str:
    mutations = []
    for i, (wt_aa, mut_aa) in enumerate(zip(wt_seq, candidate_seq)):
        if wt_aa != mut_aa: mutations.append(f"{wt_aa}{chain_id}{seq_to_pdb_map[i]}{mut_aa}")
    return ",".join(mutations) + ";" if mutations else ";"

# def compute_ddg_with_foldx(
#     seq: str, wt_seq: str, chain_id: str, seq_to_pdb_map: List[str],
#     repaired_pdb_path: str, foldx_exec: str, foldx_workdir: str, cache_dir: str,
#     number_of_runs: int, timeout_sec: int, cleanup_tmp: bool, verbose_foldx: bool
# ) -> Optional[float]:
#     import joblib
#     # MODIFIED: Passes cache_dir to get_cache_path
#     cache_path = get_cache_path(seq, cache_dir)
#     if os.path.exists(cache_path): return joblib.load(cache_path)
    
#     mutation_str = get_foldx_mutation_string(seq, wt_seq, chain_id, seq_to_pdb_map)
#     if not mutation_str or mutation_str == ";": return 0.0
    
#     tmp_dir = os.path.join(foldx_workdir, f"foldx_run_{uuid.uuid4().hex[:8]}")
#     os.makedirs(tmp_dir)
#     ddg = None
#     try:
#         repaired_base = os.path.basename(repaired_pdb_path)
#         shutil.copy2(repaired_pdb_path, os.path.join(tmp_dir, repaired_base))
#         shutil.copy2(os.path.join(foldx_workdir, "rotabase.txt"), tmp_dir)
#         with open(os.path.join(tmp_dir, "individual_list.txt"), "w") as f: f.write(mutation_str + "\n")
#         cmd = [foldx_exec, "--command=BuildModel", f"--pdb={repaired_base}", "--mutant-file=individual_list.txt", f"--numberOfRuns={number_of_runs}"]
#         proc = subprocess.run(cmd, cwd=tmp_dir, capture_output=True, text=True, timeout=timeout_sec)
#         if verbose_foldx: print(f"\n--- FOLDX Log for {seq[:10]}... ---\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}\n---------------------\n")
#         if proc.returncode == 0:
#             avg_fxout = os.path.join(tmp_dir, f"Average_{os.path.splitext(repaired_base)[0]}.fxout")
#             ddg = parse_ddg_from_fxout(avg_fxout)
#     except Exception: pass
#     finally:
#         if cleanup_tmp: shutil.rmtree(tmp_dir, ignore_errors=True)
#     if ddg is not None: joblib.dump(ddg, cache_path)
#     return ddg

# for reporting the foldx output to log

def compute_ddg_with_foldx(
    seq: str, wt_seq: str, chain_id: str, seq_to_pdb_map: List[str],
    repaired_pdb_path: str, foldx_exec: str, foldx_workdir: str, cache_dir: str,
    number_of_runs: int, timeout_sec: int, cleanup_tmp: bool, verbose_foldx: bool,
    step: int, log_file_path: str  
) -> Optional[float]:
    import joblib
    cache_path = get_cache_path(seq, cache_dir)
    if os.path.exists(cache_path): return joblib.load(cache_path)
    
    mutation_str = get_foldx_mutation_string(seq, wt_seq, chain_id, seq_to_pdb_map)
    if not mutation_str or mutation_str == ";": return 0.0
    
    tmp_dir = os.path.join(foldx_workdir, f"foldx_run_{uuid.uuid4().hex[:8]}")
    os.makedirs(tmp_dir)
    ddg = None
    try:
        repaired_base = os.path.basename(repaired_pdb_path)
        shutil.copy2(repaired_pdb_path, os.path.join(tmp_dir, repaired_base))
        shutil.copy2(os.path.join(foldx_workdir, "rotabase.txt"), tmp_dir)
        with open(os.path.join(tmp_dir, "individual_list.txt"), "w") as f: f.write(mutation_str + "\n")
        cmd = [foldx_exec, "--command=BuildModel", f"--pdb={repaired_base}", "--mutant-file=individual_list.txt", f"--numberOfRuns={number_of_runs}"]
        
        proc = subprocess.run(cmd, cwd=tmp_dir, capture_output=True, text=True, timeout=timeout_sec)

        if step == 1 and log_file_path:
            with open(log_file_path, 'a') as f:
                f.write(f"\n--- FOLDX Log for sequence {seq[:15]}... ---\n")
                f.write(f"STDOUT:\n{proc.stdout}\n")
                f.write(f"STDERR:\n{proc.stderr}\n")
                f.write("---------------------------------------------------\n")
        
        if proc.returncode == 0:
            avg_fxout = os.path.join(tmp_dir, f"Average_{os.path.splitext(repaired_base)[0]}.fxout")
            ddg = parse_ddg_from_fxout(avg_fxout)
    except Exception as e:
        if log_file_path:
            with open(log_file_path, 'a') as f:
                f.write(f"[ERROR] Exception during FoldX run for seq {seq[:15]}: {e}\n")
    finally:
        if cleanup_tmp: shutil.rmtree(tmp_dir, ignore_errors=True)
    
    if ddg is not None: joblib.dump(ddg, cache_path)
    return ddg

# class FoldXScorer(GuidedDecodingScoringFunction):
#     # MODIFIED: __init__ now accepts **kwargs to catch all extra settings
#     def __init__(self, wt_seq, chain_id, seq_to_pdb_map, repaired_pdb_path, **kwargs):
#         self.wt_seq = wt_seq
#         self.chain_id = chain_id
#         self.seq_to_pdb_map = seq_to_pdb_map
#         self.repaired_pdb_path = repaired_pdb_path
#         self.kwargs = kwargs # Store all other foldx settings
    
#     # MODIFIED: __call__ now passes the stored settings along to the compute function
#     def __call__(self, protein: ESMProtein) -> float:
#         ddg = compute_ddg_with_foldx(
#             protein.sequence, self.wt_seq, self.chain_id, self.seq_to_pdb_map, self.repaired_pdb_path, **self.kwargs
#         )
#         if ddg is None: return float("-inf")
#         score = -float(ddg)
#         return score



# class FoldXScorer(GuidedDecodingScoringFunction):
#     def __init__(self, wt_seq, chain_id, seq_to_pdb_map, repaired_pdb_path, **kwargs):
#         self.wt_seq = wt_seq
#         self.chain_id = chain_id
#         self.seq_to_pdb_map = seq_to_pdb_map
#         self.repaired_pdb_path = repaired_pdb_path
#         self.kwargs = kwargs
#         # NEW: Define the set of invalid characters
#         self.invalid_amino_acids = {'B', 'J', 'O', 'U', 'X', 'Z'}
    
#     def __call__(self, protein: ESMProtein) -> float:
#         sequence = protein.sequence

#         # MODIFIED: Expanded check for any invalid character
#         if any(char in self.invalid_amino_acids for char in sequence):
#             print(f"[SCORE] Invalid sequence with non-standard amino acid found. Score: -inf")
#             return float("-inf")

#         # The rest of the function is the same
#         ddg = compute_ddg_with_foldx(
#             sequence, self.wt_seq, self.chain_id, self.seq_to_pdb_map, self.repaired_pdb_path, **self.kwargs
#         )
#         if ddg is None: 
#             return float("-inf")
            
#         score = -float(ddg)
#         return score

# for reporting the foldx output to log

class FoldXScorer(GuidedDecodingScoringFunction):
    def __init__(self, wt_seq, chain_id, seq_to_pdb_map, repaired_pdb_path, **kwargs):
        self.wt_seq = wt_seq
        self.chain_id = chain_id
        self.seq_to_pdb_map = seq_to_pdb_map
        self.repaired_pdb_path = repaired_pdb_path
        self.kwargs = kwargs
        self.invalid_amino_acids = {'B', 'J', 'O', 'U', 'X', 'Z'}

    def __call__(self, protein: ESMProtein, step: int, log_file_path: str) -> float:
        sequence = protein.sequence

        if any(char in self.invalid_amino_acids for char in sequence):
            print(f"[SCORE] Invalid sequence with non-standard amino acid found. Score: -inf")
            return float("-inf")

        ddg = compute_ddg_with_foldx(
            sequence, self.wt_seq, self.chain_id, self.seq_to_pdb_map, 
            self.repaired_pdb_path, step=step, log_file_path=log_file_path, **self.kwargs
        )
        
        if ddg is None: 
            return float("-inf")
            
        score = -float(ddg)
        return score



def plot_ddg_history(all_scores_history: Dict, save_path: str):
    """
    Creates and saves a box plot of ddG scores over generation steps,
    with individual data points overlaid.
    """
    plot_data = []
    for step, scores in all_scores_history.items():
        for score in scores:
            if score is not None and score != float('-inf'):
                ddg = -score
                plot_data.append({'step': step, 'ddg': ddg})
    
    if not plot_data:
        print("[WARN] No valid scores to plot.")
        return

    df = pd.DataFrame(plot_data)
    
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(12, 7))
    
    sns.boxplot(x='step', y='ddg', data=df, ax=ax, palette="viridis")
    
    sns.stripplot(
        x='step', 
        y='ddg', 
        data=df, 
        ax=ax, 
        color='black', 
        alpha=0.6, 
        jitter=True  
    )
    
    ax.axhline(0, color='red', linestyle='--', linewidth=2, label='Wild-Type Stability (ΔΔG = 0)')
    ax.set_title('Distribution of Predicted ΔΔG Scores per Generation Step', fontsize=16)
    ax.set_xlabel('Generation Step', fontsize=12)
    ax.set_ylabel('Predicted ΔΔG (kcal/mol)', fontsize=12)
    ax.legend()
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\n[INFO] Plot saved to: {save_path}")
    plt.show()