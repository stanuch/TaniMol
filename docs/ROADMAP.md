# Roadmap

Current status and planned development for TaniMol.


## Phase 1 — Project setup ✅

- [x] Define research question and scope
- [x] Design project structure (directories, modules)
- [x] Write README with methodology overview
- [x] Write METHODS documentation
- [x] Set up repository (LICENSE, .gitignore)


## Phase 2 — Data acquisition

- [x] Research and select DNA repair targets from ChEMBL
- [x] Write `scripts/fetch_data.py` to download ChEMBL SQLite database dump
- [ ] Inspect the data: number of compounds per target, activity ranges, quality


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


## Phase 6 — DNA repair-specific analyses

These analyses go beyond the generic similarity pipeline and are specific to the DNA repair inhibitor domain. They require the core pipeline (Phases 2–5) to be complete first.

### Cross-target selectivity
- [ ] Identify molecules tested against multiple targets (e.g. PARP1 **and** PARP2)
- [ ] Compare activity profiles: is a potent PARP1 inhibitor also potent against PARP2?
- [ ] Determine which structural features drive selectivity between closely related targets
- [ ] Clinically relevant — olaparib vs niraparib have different selectivity profiles and side effects

### Pathway-level comparison
- [ ] Group targets by DNA repair pathway (BER: PARP1/PARP2, Checkpoint: ATR/ATM, NHEJ: DNA-PKcs)
- [ ] Compare chemical space between pathways — are inhibitors of one pathway structurally similar to another?
- [ ] High cross-pathway similarity = risk of off-target effects; low similarity = distinct binding pockets

### Synthetic lethality exploration
- [ ] Search for molecules that hit multiple targets simultaneously (dual inhibitors)
- [ ] Focus on known synthetic lethality pairs (e.g. PARP + ATR, PARP + DNA-PKcs)
- [ ] Relevant for combination therapy design in oncology

### Approved drug benchmarking
- [ ] Pull approved drugs (`max_phase = 4`) from ChEMBL for each target
- [ ] Map olaparib, niraparib, rucaparib etc. onto the Tanimoto similarity clusters
- [ ] Check whether nearest structural neighbors of approved drugs share similar activity
- [ ] Provides clinical context to cluster-level findings

### Scaffold analysis
- [ ] Extract Murcko scaffolds (RDKit) for each compound
- [ ] Count unique scaffolds per target — how chemically diverse is each target's inhibitor set?
- [ ] Identify scaffolds shared between targets (potential polypharmacology)
- [ ] Compare scaffold diversity of well-studied (PARP1) vs. underexplored (ATM) targets


## Possible extensions (not committed)

Additional ideas that might be worth exploring:

- [ ] Interactive Plotly dashboard
- [ ] Comparison of fingerprint types (how much do results change with MACCS vs Morgan?)
- [ ] Matched molecular pair analysis (systematic single-atom changes)
