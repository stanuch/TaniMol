# Changelog

All notable changes to the **TaniMol** project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.0] - 2026-04-16

### Added
- `src/visualization.py` — five new activity-focused plotting functions:
  - `plot_cluster_activity_boxplots()` — pIC50 distributions per cluster
  - `plot_activity_cliff_scatter()` — Tanimoto vs |ΔpIC50| with cliff zone overlay
  - `plot_sali_distribution()` — log-scaled SALI histogram with 95th/99th percentile markers
  - `plot_similarity_activity_density()` — hexbin density plot with Spearman ρ annotation
- `src/export.py` — results export module:
  - `export_cluster_statistics()` — per-cluster pIC50 summary to CSV
  - `export_activity_cliffs()` — cliff pairs with SMILES, similarity, SALI to CSV
  - `export_molecule_summary()` — per-molecule table with cluster ID, max SALI, cliff count to CSV
  - `export_pipeline_summary()` — full run metadata and key results to JSON
  - `export_all_figures()` — save all plots as 300 DPI PNGs
- `notebooks/07_export.ipynb` — export workflow notebook
- Centralized `COLORS` dictionary in `visualization.py` for consistent plot styling
- `RESULTS_DIR` path in `config.py`

### Changed
- `visualization.py` — unified color palette across all plots (teal/coral/amber theme)
- `visualization.py` — activity cliff scatter now uses `Rectangle` patch in data coordinates
  for accurate cliff zone alignment
- `visualization.py` — SALI histogram uses log-scaled axes and log-spaced bins to show
  the full distribution including rare high-SALI outliers
- README updated: stage badge → MVP Complete, visualization section expanded,
  project structure updated with all notebooks and export module
- ROADMAP updated: Phase 3 and Phase 4 items marked complete


## [0.7.0] - 2026-04-13

### Added
- `src/activity_analysis.py` — completed all activity analysis functions:
  - `activity_cliffs()` — vectorized detection of molecular pairs with high structural
    similarity (Tanimoto > 0.8) but dramatically different activity (|ΔpIC50| > 2.0);
    configurable thresholds, returns results sorted by most dramatic cliff first
  - `sali()` — Structure-Activity Landscape Index computation for all molecular pairs;
    returns full NxN SALI matrix and upper-triangle values vector
  - `similarity_activity_correlation()` — global SAR hypothesis test via Spearman rank
    correlation between pairwise Tanimoto similarity and |ΔpIC50|
- `notebooks/05_activity_analysis.ipynb` — completed all TODO cells:
  - Activity cliff detection with results table
  - SALI computation with summary statistics
  - Top SALI pairs extraction with SMILES inspection
  - Spearman correlation with SAR interpretation

### Changed
- `src/activity_analysis.py` — added `scipy` import for `spearmanr`


## [0.6.0] - 2026-03-29

### Added
- `src/activity_analysis.py` — within-cluster activity distribution statistics (mean, median, std, min, max per cluster)
- `notebooks/05_activity_analysis.ipynb` — activity analysis notebook (Morgan ECFP4 only)

### Fixed
- **IC50 range filter**: added `standard_value > 0 AND standard_value <= 100000` to SQL query,
  removing biologically impossible measurements (e.g. IC50 = 12M from ChEMBL data entry errors)

### Changed
- `config.py` — narrowed `TARGETS` to PARP1 only for clean single-target SAR analysis
- Renamed `notebooks/05_visualization.ipynb` to `notebooks/06_visualization.ipynb`


## [0.5.0] - 2026-03-17

### Fixed
- **Butina clustering**: replaced `squareform` with correct lower-triangle list comprehension —
  previous implementation silently produced chemically invalid clusters by mapping distances to wrong molecule pairs
- **Deduplication**: switched from keeping the most potent (first) IC50 to geometric median on
  log-scale values, eliminating best-case bias in the activity dataset
- **Fingerprint pipeline**: invalid SMILES are now filtered before fingerprint generation,
  preserving DataFrame alignment and preventing silent `None` propagation downstream
- **Stereochemistry**: enforced `isomericSmiles=True` in SMILES canonicalization to prevent
  collapse of enantiomers with distinct biological activity
- **ChEMBL query**: restricted `standard_relation` to `('=', '~')`, excluding censored
  measurements (`>`, `<`) that were previously treated as exact values
- **Activity types**: restricted `ACTIVITY_TYPES` to `IC50` only, removing naive pooling
  with `Ki` (a distinct thermodynamic measure)

### Changed
- Similarity heatmap now uses UPGMA (average linkage) instead of Ward linkage —
  Ward assumes Euclidean geometry which is invalid for Tanimoto distances on binary fingerprints
- `standardize_molecules` refactored from `iterrows()` to `pandas.apply()` for performance

### Added
- `visualization.py` — plotting utilities: similarity distribution, cluster size distribution,
  top-N cluster heatmap, and full similarity heatmap (sanity check)
- `notebooks/05_visualization.ipynb` — visual analysis of clustering results

### Chore
- Pinned all dependency versions in `environment.yml`


## [0.4.0] - 2026-03-16

### Added
- `src/fingerprints.py` — fingerprint generation module:
  - `mol_from_smiles()` — parse SMILES strings safely
  - `generate_morgan_fp()`, `generate_maccs_fp()`, `generate_rdkit_fp()` — generators for 3 fingerprint types
  - `add_fingerprints()` — batch generation for datasets
  - `save_fingerprints()`, `load_fingerprints()` — persist arrays to `.npy`
- `src/similarity.py` — similarity matrix computation:
  - `calculate_similarity_matrix()` — optimized underlying NumPy `np.dot` implementation for fast Tanimoto computation
- `src/visualization.py` — plotting functions:
  - `plot_similarity_heatmap()` — generates a Tanimoto heatmap sorted by UPGMA hierarchical clustering
- `notebooks/02_fingerprints.ipynb` — workflow for generating and storing dataset fingerprints
- `notebooks/03_similarity.ipynb` — workflow for calculating matrices and performing EDA/sanity checks

### Changed
- `src/config.py` — added fingerprint storage paths (`FINGERPRINTS_DIR`, `MORGAN_FP_PATH`, etc.)


## [0.3.0] - 2026-03-13

### Added
- `src/preprocessing.py` — full preprocessing pipeline with functions:
  - `fetch_activity_data()` — query ChEMBL for activity data across all configured targets
  - `drop_missing_values()` — remove rows with missing `standard_value`
  - `validate_smiles()` — remove rows with invalid SMILES (RDKit validation)
  - `standardize_molecules()` — strip salts, neutralize charges, canonicalize tautomers
  - `deduplicate()` — keep best (lowest) IC50 per (target, SMILES) pair
  - `compute_pic50()` — fill missing pchembl_value using −log₁₀(IC50 × 10⁻⁹)
  - `save_cleaned_data()` — export to CSV with per-target summary
- `notebooks/01_analysis.ipynb` — analysis notebook running the preprocessing pipeline
- `tqdm` added to project dependencies

### Changed
- `src/config.py` — added output file path configuration (`OUTPUT_FILE_NAME`, `OUTPUT_PATH`)
- `pyproject.toml` — added `tqdm` to dependencies, fixed `build-backend` to `setuptools.build_meta`


## [0.2.0] - 2026-02-19

### Added
- `scripts/fetch_data.py` — downloads the full ChEMBL SQLite database dump from the EBI FTP server, extracts the `.db` file into `data/raw/`, and cleans up temporary files. Includes download progress bar and skip-if-exists logic.

### Changed
- Updated README "Data acquisition" section to reflect the SQLite dump approach instead of the ChEMBL API
- Added ChEMBL database files (`*.db`, `*.tar.gz`) to `.gitignore`


## [0.1.0] - 2026-02-18

### Added
- Project structure: `src/`, `scripts/`, `notebooks/`, `data/`, `tests/`, `results/`, `docs/`
- README with research context, methodology overview, pipeline diagram, and project structure
- METHODS.md with detailed documentation of all computational methods
- ROADMAP.md with development phases and planned features
- USAGE.md with step-by-step instructions and troubleshooting
- Project logo (`docs/img/tanimol_logo.png`)
- MIT License
- `.gitignore` for Python projects