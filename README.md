# AlphaFold3-Guided Structural Modeling and Molecular Docking of Marine Sulfated Compounds

## MSc Bioinformatics & Data Science | Project

### Overview
Comparative computational study of two marine sulfated compound classes —
ulvan oligosaccharides (green algae) and sulfated sterols (marine sponges) —
using AlphaFold3 structure prediction, molecular docking (AutoDock Vina),
and molecular dynamics simulation (GROMACS), against diabetes (DPP-4) and
Alzheimer's disease (BACE-1) targets.

### Status: Phase 0 — Setup Complete

### Environment
- Conda env: `marine_meta`
- Tools: AutoDock Vina 1.1.2, Open Babel 3.1.0, RDKit 2024.09.3, PyMOL

### Compounds
| Compound | Source | Status |
|---|---|---|
| Ulvan oligosaccharides | *Ulva* sp. (green algae) | Pending |
| Sulfated sterols | Marine sponges | Pending |

### Targets
| Target | Disease | PDB ID |
|---|---|---|
| DPP-4 | Type 2 Diabetes | TBD - verify |
| BACE-1 | Alzheimer's | TBD - verify |

### Pipeline
AlphaFold3 -> AutoDock Vina (cross-validation) -> GROMACS MD -> ADMET

### Author
Ganapathirajan — MSc Bioinformatics & Data Science
