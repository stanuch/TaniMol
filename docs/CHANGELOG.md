# Changelog

All notable changes to the **TaniMol** project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


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