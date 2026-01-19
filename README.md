<p align="center">
  <img src="docs/img/tanimol_logo.png" width="600" title="TaniMol Logo" alt="TaniMol Logo">
</p>


**TaniMol** will be an open-source Python tool designed for chemical space exploration and molecular clustering based on Tanimoto similarity metrics. Built upon the RDKit framework, it will automate the generation of fingerprints and similarity matrices to identify structural relationships within large chemical libraries.

It will aim to provide a reproducible and efficient pipeline for selecting representative compounds, facilitating high-throughput analysis in early-stage drug discovery.

## Theoretical Workflow

1.  **Data Ingestion & Sanitization**
    The system will accept chemical structures (SMILES, SDF) and strictly validate them. It will handle salt stripping, tautomer standardization, and duplicate removal to ensure data integrity before any calculation begins.

2.  **Fingerprint Generation**
    TaniMol will map chemical structures into high-dimensional bit vectors (Morgan Fingerprints/ECFP4). This transforms chemical intuition into a mathematical format suitable for vector operations.

3.  **Similarity Matrix Calculation**
    Using the **Tanimoto Coefficient** (Jaccard Index), the tool will compute an *N x N* similarity matrix. This step allows the quantification of "likeness" between any two molecules in the set.

4.  **Unsupervised Clustering**
    Based on the distance matrix ($1 - Tanimoto$), molecules will be grouped into clusters. The primary algorithm will be **Butina clustering**, optimized for chemical datasets, with optional support for K-Means and Hierarchical clustering.

5.  **Visualization & Reporting**
    Finally, the high-dimensional data will be projected into 2D space using dimensionality reduction techniques (t-SNE or UMAP), allowing users to visually inspect "islands" of chemical activity.

## Methodological Approach

This project implements a standard chemoinformatics pipeline focused on unsupervised learning. The methodology is built upon three core pillars chosen for their proven effectiveness in drug discovery campaigns:

### 1. Molecular Representation: Morgan Fingerprints (ECFP)
Chemical structures are vectorized using **Morgan Fingerprints** (equivalent to Extended-Connectivity Fingerprints, ECFP4/ECFP6) with a radius of 2 or 3 and a bit-length of 1024 or 2048.

### 2. Similarity Metric: Tanimoto Coefficient
To quantify the relationship between binary fingerprint vectors, the **Tanimoto Coefficient** ($T_c$) is utilized:

$$T(A, B) = \frac{c}{a + b - c}$$

Where $c$ is the count of common set bits, and $a$ and $b$ are the set bits in molecule $A$ and $B$, respectively.

### 3. Clustering Strategy: Butina Algorithm
For grouping compounds, the project prioritizes the **Butina clustering algorithm** (sphere exclusion) over general-purpose methods like K-Means.

## Planned Features

* **Robust I/O:** Support for .csv, .smi, and .sdf files with error logging for malformed structures.
* **Configurable Fingerprints:** Support for Morgan (ECFP), MACCS Keys, and Topological Torsion fingerprints.
* **Efficient Matrix Operations:** Utilizing NumPy and RDKit's bulk similarity functions to handle datasets of 10k+ molecules efficiently.
* **Interactive Visualization:** Integration with Plotly to generate interactive scatter plots of the chemical space.
* **Scaffold Analysis:** Automatic extraction of Murcko Scaffolds for each generated cluster to identify the core substructures.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.