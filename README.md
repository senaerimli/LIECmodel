# LIECmodel

Genome-scale metabolic model of the large intestine epithelial cell (LIEC) for ulcerative colitis (UC):
consensus Healthy and UC models, transcriptomic and constraint-based drug-response analysis.

> Status: work in progress. Associated manuscript: (preprint/DOI to be added).

## What is in this repository
- `Models/`: final models (SBML `.xml` and small `.mat`)
- `MainCodes/`: scripts in the order given by their part numbers (e.g., `part01_...`, `part07_...`)
- `Functions/`: reusable functions
- `Tasks/`: metabolic task panel (272 tasks) with source credits
- `Data/`: small curated inputs and GEO download scripts (raw GEO data is not stored here)
- `Results/`: small result tables; large files are on Zenodo (link below)
- `docs/`: MEMOTE report, methods notes

## Requirements
See `environment.md`. Test with a free solver and state which results depend on a commercial solver.

## How to reproduce
1. Install MATLAB and the COBRA Toolbox (version in `environment.md`); run `initCobraToolbox`.
2. Clone this repository.
3. Run the scripts in `MainCodes/` in part-number order.

## Data and large files
Large models and results: Zenodo DOI (to be added).

## Citation
See `CITATION.cff` (to be added at first release).

## License
Code: (choose MIT or Apache-2.0). Models and data: CC BY 4.0 unless a source requires otherwise.

## Credits
Harvey/Harvetta whole-body models and the OrganCompendium; metabolic task list sources (list references).
