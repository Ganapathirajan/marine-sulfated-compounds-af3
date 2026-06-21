"""
Build SMILES for ulvan disaccharide repeat units: A3S, B3S, U3S
Run this in your `marine_meta` conda env on WSL.

Structures (confirmed from literature):
  A3S: beta-D-GlcA-(1->4)-alpha-L-Rha3S
  B3S: alpha-L-IdoA-(1->4)-alpha-L-Rha3S
  U3S: beta-D-Xyl-(1->4)-alpha-L-Rha3S

NOTE: These SMILES are manually constructed from published monosaccharide
linkage descriptions. They are NOT pulled from a verified database entry.
Before using in docking, visually compare the 2D RDKit-rendered structure
against the structure figures in the source papers (JACS 2025, MDPI 2024
Marine Drugs paper, Carbohydrate Research 1998 Lahaye paper) to confirm
stereochemistry and linkage position are correct. RDKit sanitization only
confirms the SMILES is chemically VALID -- not that it matches literature.
"""

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors

# --- Disaccharide SMILES (manually built, alpha-L-Rha 3-O-sulfate as aglycone unit) ---
structures = {
    "A3S": {
        "smiles": "O=S(=O)(O)O[C@@H]1[C@@H](O)[C@H](O)[C@@H](C)O[C@@H]1O[C@@H]1O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]1O",
        "desc": "beta-D-GlcA-(1->4)-alpha-L-Rha3S"
    },
    "B3S": {
        "smiles": "O=S(=O)(O)O[C@@H]1[C@@H](O)[C@H](O)[C@@H](C)O[C@@H]1O[C@H]1O[C@@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]1O",
        "desc": "alpha-L-IdoA-(1->4)-alpha-L-Rha3S"
    },
    "U3S": {
        "smiles": "O=S(=O)(O)O[C@@H]1[C@@H](O)[C@H](O)[C@@H](C)O[C@@H]1O[C@@H]1OC[C@H](O)[C@@H](O)[C@H]1O",
        "desc": "beta-D-Xyl-(1->4)-alpha-L-Rha3S"
    },
}

print(f"{'Unit':<6} {'Valid':<7} {'MW':<10} {'Formula':<15} Description")
print("-" * 70)

for name, data in structures.items():
    mol = Chem.MolFromSmiles(data["smiles"])
    if mol is None:
        print(f"{name:<6} {'FAIL':<7} -- could not parse SMILES, needs fixing")
        continue
    mw = Descriptors.MolWt(mol)
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    print(f"{name:<6} {'OK':<7} {mw:<10.2f} {formula:<15} {data['desc']}")

    # Save 2D depiction for manual visual check against literature figures
    AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol, f"compounds/ulvan/{name}_2D.png", size=(500, 500))

    # Save as SDF (3D) for docking pipeline
    mol3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol3d, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol3d)
    writer = Chem.SDWriter(f"compounds/ulvan/{name}_3d.sdf")
    writer.write(mol3d)
    writer.close()

    # Save canonical SMILES to text file
    with open(f"compounds/ulvan/{name}.smi", "w") as f:
        f.write(Chem.MolToSmiles(mol) + f"\t{name}\n")

print("\nDone. Check compounds/ulvan/ for PNG (visual check), SDF (3D), SMI (SMILES) files.")
print("\n*** MANDATORY: open each PNG and compare against published structure figures ***")
print("*** before trusting these for docking. Stereochemistry on manually-built SMILES ***")
print("*** is the most likely source of error and MUST be visually verified. ***")
