"""
Clean target PDB files for docking:
- Extract protein-only structure (no water, no inhibitor, no sugars/ions)
- Extract reference inhibitor separately (used to define binding pocket
  grid box center for docking)

Run in marine_meta conda env on WSL, from project root.
"""

def clean_pdb(input_path, protein_out, ligand_out, ligand_code):
    protein_lines = []
    ligand_lines = []
    with open(input_path) as f:
        for line in f:
            if line.startswith(("ATOM", "TER")):
                protein_lines.append(line)
            elif line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname == ligand_code:
                    ligand_lines.append(line)
            elif line.startswith(("HEADER", "TITLE", "REMARK", "END")):
                pass  # skip metadata for clean coordinate files

    with open(protein_out, "w") as f:
        f.writelines(protein_lines)
        f.write("END\n")

    with open(ligand_out, "w") as f:
        f.writelines(ligand_lines)
        f.write("END\n")

    print(f"{input_path}: {len(protein_lines)} protein atom lines, "
          f"{len(ligand_lines)} ligand ({ligand_code}) atom lines")

# DPP-4: reference inhibitor code = 715
clean_pdb(
    "targets/dpp4/1X70.pdb",
    "targets/dpp4/1X70_protein_clean.pdb",
    "targets/dpp4/1X70_ligand_715.pdb",
    "715"
)

# BACE-1: reference inhibitor code = 3BN
clean_pdb(
    "targets/bace1/2B8V.pdb",
    "targets/bace1/2B8V_protein_clean.pdb",
    "targets/bace1/2B8V_ligand_3BN.pdb",
    "3BN"
)

print("\nDone. Check targets/dpp4/ and targets/bace1/ for cleaned files.")
print("NOTE: 1X70 has 2 chains (A, B) -- both kept here. Decide later if")
print("you dock against the full dimer or a single chain (common practice")
print("is single chain for simpler docking; check literature precedent).")
