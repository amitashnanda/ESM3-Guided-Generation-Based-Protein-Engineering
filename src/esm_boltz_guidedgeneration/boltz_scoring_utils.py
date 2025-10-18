"""This script implements a subprocess to run boltz predict for protein-ligand
affinity prediction given a protein sequence and ligand SMILES string.

See example outputs in boltz_outputs/.

Author: Anna Su <anna.su@yale.edu>
Date: 2025-10-17
"""

import os
import json
import subprocess
import shutil
from datetime import datetime
from pathlib import Path


def get_boltz_affinity(seq, smiles, out_dir):
     
    # Create a yaml file with the seq and smiles
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    yaml_path = os.path.join(out_dir, f"protein_ligand.yaml")
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

    # Run boltz predict
    boltz_exe = shutil.which("boltz")
    cmd = [boltz_exe, "predict", yaml_path, "--out_dir", out_dir, "--use_msa_server"]
    _ = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Get the affinity prediction from the output json
    json_files = list(Path(out_dir).rglob("affinity_*.json"))
    json_path = max(json_files, key=lambda p: p.stat().st_mtime)
    with open(json_path, "r") as jf:
        data = json.load(jf)
    print(f"Boltz predicted affinity: {data['affinity_pred_value']}")
    return float(data["affinity_pred_value"])

# Example usage
# seq = "ACNYTCGSNVYSSSQVDAYLATGYKLHEDGVTVGSNSYPHKYNNMEGFQFRVSSGYYEWPILDSGRTYSGGSPGADRVVFNENGILAGVILHYGASGNNFVECT"
# smiles = "NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O"
# out_dir = f"boltz_outputs/run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
# a = get_boltz_affinity(seq, smiles, out_dir)