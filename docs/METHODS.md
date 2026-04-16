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


## Clustering

Molecules are grouped using **Butina clustering** (sphere-exclusion). The algorithm ranks molecules by neighbor count and assigns each molecule to a cluster if it falls within a Tanimoto distance threshold (default: 0.4, i.e. similarity ≥ 0.6) of the cluster centroid. Molecules with no neighbors within the threshold remain as singletons.

UPGMA (average linkage) hierarchical clustering is used exclusively for heatmap ordering — it sorts the similarity matrix to surface diagonal "islands" of structural analogs. Ward linkage was intentionally avoided because it assumes Euclidean geometry, which is invalid for Tanimoto distances on binary fingerprints.


## Activity Analysis

**Activity cliffs** are defined as pairs of molecules with Tanimoto similarity > 0.8 and |ΔpIC50| > 2.0 (corresponding to >100-fold potency difference). Detection is fully vectorized via NumPy boolean masking on the precomputed similarity and |ΔpIC50| matrices.

**SALI (Structure-Activity Landscape Index)** provides a continuous score for every molecular pair:

$$SALI(i,j) = \frac{|\Delta pIC50|}{1 - T(i,j)}$$

SALI amplifies cases where structurally very similar molecules have very different activities. Division by zero (identical fingerprints, T = 1.0) is handled by setting SALI to 0. Note: Morgan ECFP4 does not encode chirality by default, so enantiomers produce T = 1.0 despite potentially significant activity differences — a known limitation documented in the project findings.

**Spearman rank correlation** between pairwise Tanimoto similarity and |ΔpIC50| (upper triangle only, to avoid counting pairs twice) serves as a global SAR hypothesis test. A negative ρ confirms that structurally similar molecules tend to have similar biological activity.
