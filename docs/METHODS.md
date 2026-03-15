# Methods

Detailed description of the computational methods used in TaniMol. This document supplements the overview in [README.md](../README.md) with implementation details, parameter choices, and references.


## Data Source and Filtering

Bioactivity data comes from the full ChEMBL SQLite dump (version 36) to bypass API pagination limits for the ~11,000 human DNA repair protein inhibitors to be processed. 

Target selection required involvement in BER, HR, NHEJ, MMR, or related pathways, with at least 50 compounds having IC50 data. I've prioritized clinically relevant targets like PARP1, PARP2, ATR, ATM, and DNA-PKcs.

I filter the raw assay data strictly for Binding (B) or Functional (F) assays with a ChEMBL confidence score ≥ 7. I select only exact (`=`) IC50 or Ki values measured in nM for single-protein targets.


## Molecular Preprocessing

SMILES standardization relies on RDKit's `MolStandardize` module. First, I parse the string, strip salts and non-covalent fragments, neutralize formal charges where reasonable, and iterate tautomers to pick the canonical limit. Duplicate entries (identical canonical SMILES for the same target) are resolved by keeping the measurement with the lowest IC50. 

All IC50 values (nM) are converted to pIC50 ($-\log_{10}(IC50 \times 10^{-9})$) before any statistical models are applied.


## Molecular Fingerprints

Primary focus is on Morgan fingerprints (ECFP). I settled on radius 2 (ECFP4 standard) hashed to 2048 bits as the default. This sparse representation keeps the median Tanimoto similarity across the ~11,000 compounds tightly around 0.13.

For baseline comparison, I also generate 166-bit MACCS Keys. Due to their low dimensionality and rigid rule-based nature, MACCS keys yield heavily inflated similarity distributions in the dataset (median similarity jumps to 0.52 compared to Morgan's 0.13, with ~60% of all pairs scoring above 0.5). 2048-bit RDKit topological path fingerprints are included as a middle-ground benchmark.


## Similarity Calculation

Pairwise Tanimoto similarity was calculated using the formula $T = \frac{c}{a+b-c}$ between fingerprint vectors. 

Given the scale of my data ($11,000 \times 11,000$ matrix yields over 121 million pairs), evaluating distances row-by-row with `scikit-learn`'s `pairwise_distances` or `scipy.spatial.distance.pdist` proved too slow for rapid iteration. I used vectorized NumPy dot-product approach. This brought computation time for the entire matrix down to roughly 0.5 seconds on a standard desktop CPU (from ~1:47 minutes).

For algorithms requiring actual distances, I apply $D = \text{clip}(1 - T, 0, 1)$ to prevent floating-point underflow.
