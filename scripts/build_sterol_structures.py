"""
Build 3D structures + 2D visual check for sulfated sterols.
FINAL compound set (Contignasterol replaced -- it had zero sulfate groups,
didn't fit the "sulfate functional group" comparative theme).

Halistanol sulfate : PubChem CID 73361 (neutral, 3 sulfates)
Solomonsterol A     : PubChem CID 91932303 (anion form) -- manually
                       protonated here ([O-] -> OH on sulfates) for docking,
                       since neutral/protonated form needed for AutoDock Vina.
                       Verify protonation state assumption against your
                       docking protocol's pH setting before final use.

Run in marine_meta conda env on WSL.
"""

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors

compounds = {
    "halistanol_sulfate": {
        "smiles": "C[C@H](CC[C@H](C)C(C)(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2C[C@@H]([C@@H]4[C@@]3(C[C@@H]([C@H](C4)OS(=O)(=O)O)OS(=O)(=O)O)C)OS(=O)(=O)O)C",
        "cid": 73361,
        "source": "PubChem CID 73361 (neutral form, as deposited)"
    },
    "solomonsterol_A": {
        "smiles": "C[C@H](CCCOS(=O)(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(C[C@@H]([C@H](C4)OS(=O)(=O)O)OS(=O)(=O)O)C)C",
        "cid": "91932303 (anion, manually protonated for docking)",
        "source": "PubChem CID 91932303, sulfate groups protonated O- -> OH"
    },
}

print(f"{'Compound':<20} {'Valid':<7} {'MW':<10} {'Formula':<15} Source")
print("-" * 80)

for name, data in compounds.items():
    mol = Chem.MolFromSmiles(data["smiles"])
    if mol is None:
        print(f"{name:<20} FAIL -- could not parse")
        continue
    mw = Descriptors.MolWt(mol)
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    print(f"{name:<20} {'OK':<7} {mw:<10.2f} {formula:<15} {data['source']}")

    AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol, f"compounds/sulfated_sterols/{name}_2D.png", size=(600, 600))

    mol3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol3d, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol3d)
    writer = Chem.SDWriter(f"compounds/sulfated_sterols/{name}_3d.sdf")
    writer.write(mol3d)
    writer.close()

    with open(f"compounds/sulfated_sterols/{name}.smi", "w") as f:
        f.write(Chem.MolToSmiles(mol) + f"\t{name}\n")

print("\nDone. PNG/SDF/SMI saved in compounds/sulfated_sterols/")
print("Visual check: both should show 4-ring steroid core + 3 sulfate groups (-OSO3H) each")
print("Expected MW: Halistanol sulfate ~688.9 | Solomonsterol A (protonated) should be close to 615.8 + 3H ~ 618.8")
