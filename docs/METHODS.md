# Methods

Detailed description of the computational methods used in TaniMol. This document supplements the overview in [README.md](../README.md) with implementation details, parameter choices, and references.


## Data Source

Bioactivity data is obtained from the [ChEMBL database](https://www.ebi.ac.uk/chembl/) (version 34 or later) via the `chembl_webresource_client` Python API.

### Target selection

The project focuses on inhibitors of DNA repair proteins. Target selection criteria:
- Protein must be involved in a DNA repair pathway (BER, HR, NHEJ, MMR, or related)
- ChEMBL must contain at least ~50 compounds with IC50 data for the target
- Preference for targets with clinical relevance (approved or clinical-stage inhibitors exist)

### Data filtering

Raw ChEMBL bioactivity data is filtered before use:

| Filter | Value | Reason |
|--------|-------|--------|
| Assay type | Binding (B) or Functional (F) | Exclude ADMET and unclassified assays |
| Confidence score | ≥ 7 | Only high-confidence target assignments |
| Standard type | IC50 or Ki | Most common and comparable potency measures |
| Standard units | nM | Ensures unit consistency |
| Standard relation | `=` | Exclude inequalities (>, <, ~) |
| Target type | Single protein | Exclude cell-based or multi-target assays |


## Molecular Preprocessing

### SMILES standardization

Molecules are processed using RDKit:

1. **Parsing** - SMILES strings are parsed into RDKit `Mol` objects. Molecules that fail parsing are logged and discarded.
2. **Salt stripping** - Using `rdkit.Chem.SaltRemover`, counter-ions, solvents, and other non-covalent fragments are removed. Only the largest fragment is retained.
3. **Charge neutralization** - Formal charges are neutralized where chemically reasonable (e.g. carboxylate anions → carboxylic acids).
4. **Tautomer canonicalization** - Using `rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator`, each molecule is converted to its canonical tautomeric form.
5. **Canonical SMILES** - The cleaned molecule is written back to SMILES using `Chem.MolToSmiles(mol, canonical=True)`.

### Deduplication

After standardization, duplicate molecules (identical canonical SMILES for the same target) are detected. When duplicates exist, the entry with the **lowest IC50** (highest potency) is retained, as it represents the most optimistic measurement.

### Activity conversion

IC50 values (in nM) are converted to **pIC50**:

$$pIC50 = -\log_{10}(IC50 \times 10^{-9})$$

This transformation has two advantages:
- Higher values correspond to more potent compounds (more intuitive)
- The logarithmic scale compresses the range, making statistical comparisons more meaningful (IC50 values can span from 0.1 nM to 100,000 nM)


## Molecular Fingerprints

### Morgan fingerprints (ECFP)

The primary fingerprint type. Each atom in a molecule is assigned an initial identifier based on its properties (element, charge, degree, etc.). These identifiers are then iteratively updated by incorporating information from neighboring atoms up to a specified radius.

**Default parameters:**
- Radius: 2 (equivalent to ECFP4, where 4 = diameter = 2 × radius)
- Number of bits: 2048
- Features: standard Morgan (atom invariants based on atomic number, degree, etc.)

**Why radius 2:** This captures enough substructural context for meaningful similarity while avoiding the fingerprint "saturation" that occurs at higher radii (where too many bits are set and all molecules look similar).

**Why 2048 bits:** A balance between resolution (more bits = fewer hash collisions = more accurate similarity) and memory/speed. 1024 bits is also commonly used, but 2048 is the more standard choice in recent literature.

### Alternative fingerprints

For comparison, the project also supports:

- **MACCS Keys** (166 bits) - Predefined structural patterns (e.g. "contains a 6-membered ring", "contains a nitrogen atom bonded to two carbons"). Less granular than Morgan but based on expert-curated features.
- **RDKit topological fingerprints** (2048 bits) - Path-based fingerprints that enumerate all linear paths of atoms up to a given length.

Comparing results across fingerprint types can reveal how sensitive the analysis is to the choice of molecular representation.


## Similarity Calculation

### Tanimoto coefficient

For two binary fingerprints A and B:

$$T(A, B) = \frac{|A \cap B|}{|A \cup B|} = \frac{c}{a + b - c}$$

where:
- *a* = number of bits set in fingerprint A
- *b* = number of bits set in fingerprint B
- *c* = number of bits set in both A and B

Properties:
- Range: [0, 1]
- T = 1.0 when fingerprints are identical
- T = 0.0 when fingerprints share no bits
- Symmetric: T(A, B) = T(B, A)
- Not a proper distance metric (does not satisfy triangle inequality), but 1−T is commonly used as a distance

### Distance matrix

The distance matrix used for clustering is computed as:

$$D(A, B) = 1 - T(A, B)$$

This is not a true metric distance (it violates the triangle inequality), but it works well in practice for Butina clustering and hierarchical clustering with average linkage.


## Clustering

### Butina clustering

The primary clustering method. Algorithm:

1. Compute the number of neighbors for each molecule (neighbors = molecules within distance threshold)
2. Sort molecules by neighbor count (descending)
3. The molecule with the most neighbors becomes a **centroid**; all its neighbors join the cluster
4. Remove assigned molecules and repeat
5. Molecules with no neighbors within the threshold become **singletons**

**Default distance threshold:** 0.4 (equivalent to Tanimoto similarity ≥ 0.6). This is a commonly used cutoff in chemoinformatics — molecules with Tanimoto ≥ 0.6 using ECFP4 fingerprints are generally considered to share a common scaffold.

The threshold is the single most important parameter. Lower values produce more, smaller clusters; higher values produce fewer, larger clusters.

### Hierarchical clustering

Applied as a secondary method using scipy's `linkage` function:
- **Linkage method:** Ward (minimizes within-cluster variance) or average (UPGMA)
- **Input:** condensed distance matrix from Tanimoto distances
- **Output:** dendrogram + flat cluster assignments at a chosen cut height

Hierarchical clustering produces a dendrogram which is useful for visualizing the overall relationship structure, even when Butina is used for the actual cluster assignments.


## Activity Analysis

### Activity cliff detection

An activity cliff is a pair of molecules that are structurally similar but have significantly different biological activity. Formally:

A pair (i, j) is an activity cliff if:
- T(i, j) ≥ similarity threshold (default: 0.8)
- |pIC50(i) − pIC50(j)| ≥ activity threshold (default: 2.0, i.e. 100-fold potency difference)

### SALI (Structure-Activity Landscape Index)

For each pair of molecules (i, j):

$$SALI(i, j) = \frac{|pIC50_i - pIC50_j|}{1 - T(i, j)}$$

SALI amplifies cases where very similar molecules have very different activity. A pair with T = 0.95 and ΔpIC50 = 3.0 gets SALI = 60, while a pair with T = 0.5 and ΔpIC50 = 3.0 gets SALI = 6. This makes SALI useful for ranking activity cliffs by severity.

Note: SALI is undefined when T = 1.0 (identical fingerprints). These pairs are excluded.

**Reference:** Guha, R.; Van Drie, J. H. *J. Chem. Inf. Model.* 2008, 48, 646–658.

### Similarity-activity correlation

The overall relationship between structural similarity and activity similarity is quantified using:

- **Spearman rank correlation** between pairwise Tanimoto values and pairwise |ΔpIC50| values
- Expected result: negative correlation (higher similarity → smaller activity difference)
- Statistical significance assessed via p-value


## Visualization

### Chemical space projection

High-dimensional fingerprint vectors are projected to 2D using:

- **t-SNE** (t-distributed Stochastic Neighbor Embedding) - preserves local structure. Good for revealing clusters but distances between distant clusters are not meaningful.
- **UMAP** (Uniform Manifold Approximation and Projection) - preserves both local and global structure better than t-SNE. Faster on large datasets.

Points are colored by target protein or by pIC50 value to overlay biological information onto the structural map.

### Other plots

| Plot | Purpose |
|------|---------|
| Similarity heatmap | Full N×N Tanimoto matrix as a color-coded grid |
| Cluster boxplots | pIC50 distribution per cluster (box or violin plot) |
| Activity cliff scatter | Tanimoto vs. \|ΔpIC50\| for all pairs, cliffs highlighted |
| SALI network | Graph where nodes = molecules, edges = high-SALI pairs |


## References

- Bajusz, D.; Rácz, A.; Héberger, K. Why is Tanimoto index an appropriate choice for fingerprint-based similarity calculations? *J. Cheminf.* 2015, 7, 20.
- Butina, D. Unsupervised Data Base Clustering Based on Daylight's Fingerprint and Tanimoto Similarity. *J. Chem. Inf. Comput. Sci.* 1999, 39, 747–750.
- Guha, R.; Van Drie, J. H. Structure-Activity Landscape Index: Identifying and Quantifying Activity Cliffs. *J. Chem. Inf. Model.* 2008, 48, 646–658.
- Rogers, D.; Hahn, M. Extended-Connectivity Fingerprints. *J. Chem. Inf. Model.* 2010, 50, 742–754.
