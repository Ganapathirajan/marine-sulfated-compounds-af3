"""
Extract protein sequence (FASTA) from cleaned PDB files.
Needed as input for AlphaFold Server (alphafoldserver.com).

Requires biopython. Install if missing:
  pip install biopython

Run in marine_meta conda env on WSL.
"""

from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder

def extract_fasta(pdb_path, chain_id, out_fasta, seq_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", pdb_path)
    ppb = PPBuilder()

    seq = ""
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for pp in ppb.build_peptides(chain):
                    seq += str(pp.get_sequence())

    if not seq:
        print(f"WARNING: no sequence extracted for chain {chain_id} in {pdb_path}")
        return

    with open(out_fasta, "w") as f:
        f.write(f">{seq_id}\n{seq}\n")

    print(f"{seq_id}: {len(seq)} residues -> {out_fasta}")
    print(seq)
    print()

# DPP-4 (1X70) - using chain A (decide: single chain for docking simplicity)
extract_fasta("targets/dpp4/1X70_protein_clean.pdb", "A",
              "targets/dpp4/1X70_chainA.fasta", "DPP4_1X70_chainA")

# BACE-1 (2B8V) - single chain A
extract_fasta("targets/bace1/2B8V_protein_clean.pdb", "A",
              "targets/bace1/2B8V_chainA.fasta", "BACE1_2B8V_chainA")

print("Done. Copy the sequence(s) above into AlphaFold Server (alphafoldserver.com)")
print("as the protein input, paired with each ligand SMILES, for each compound-target run.")
