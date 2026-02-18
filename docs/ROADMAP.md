# Roadmap

Current status and planned development for TaniMol.


## Phase 1 — Project setup ✅

- [x] Define research question and scope
- [x] Design project structure (directories, modules)
- [x] Write README with methodology overview
- [x] Write METHODS documentation
- [x] Set up repository (LICENSE, .gitignore)


## Phase 2 — Data acquisition

- [ ] Research and select DNA repair targets from ChEMBL
- [ ] Write `scripts/fetch_data.py` to pull bioactivity data via ChEMBL API
- [ ] Download and save per-target CSVs to `data/raw/`
- [ ] Inspect the data: number of compounds per target, activity ranges, quality
- [ ] Include a fallback CSV in the repo for offline use


## Phase 3 — Core pipeline

- [ ] `src/preprocessing.py` — SMILES validation, salt stripping, standardization, pIC50 conversion
- [ ] `src/fingerprints.py` — Morgan (ECFP4), MACCS, RDKit fingerprint generation
- [ ] `src/similarity.py` — Tanimoto similarity matrix + distance matrix
- [ ] `src/clustering.py` — Butina clustering + hierarchical clustering
- [ ] `src/activity_analysis.py` — Activity cliff detection, SALI, correlation statistics
- [ ] `src/visualization.py` — Heatmaps, chemical space maps, boxplots, cliff plots, SALI network


## Phase 4 — Analysis notebook

- [ ] `notebooks/01_analysis.ipynb` — Full end-to-end analysis
  - Data loading and EDA
  - Chemical space visualization (t-SNE/UMAP)
  - Similarity matrix and clustering
  - SAR analysis and activity cliffs
  - Conclusions


## Phase 5 — Testing and polish

- [ ] Unit tests for preprocessing, fingerprints, and similarity modules
- [ ] Run the full pipeline end-to-end, verify all outputs
- [ ] Export key figures to `results/`
- [ ] Final README review


## Possible extensions (not committed)

These are ideas that might be worth exploring if the core project works well:

- [ ] Multi-target selectivity analysis (compounds active against >1 target)
- [ ] Scaffold decomposition (Murcko scaffolds per cluster)
- [ ] Interactive Plotly dashboard
- [ ] Comparison of fingerprint types (how much do results change with MACCS vs Morgan?)
- [ ] Matched molecular pair analysis (systematic single-atom changes)
