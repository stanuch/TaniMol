# DNA Repair Protein Targets

Selected targets for the TaniMol analysis, all from *Homo sapiens*.


## Target list

| ChEMBL ID | Protein | Short name | Repair pathway | Records (IC50/Ki, ≥7 conf.) |
|---|---|---|---|---:|
| CHEMBL3105 | Poly [ADP-ribose] polymerase 1 | PARP1 | BER | 5,898 |
| CHEMBL5366 | Poly [ADP-ribose] polymerase 2 | PARP2 | BER | 1,126 |
| CHEMBL5024 | Serine/threonine-protein kinase ATR | ATR | Checkpoint | 5,110 |
| CHEMBL3797 | Serine-protein kinase ATM | ATM | Checkpoint | 1,672 |
| CHEMBL3142 | DNA-dependent protein kinase catalytic subunit | DNA-PKcs | NHEJ | 3,078 |


## Why these targets

These five proteins were chosen because they span three distinct DNA repair pathways while being part of the same overall DNA damage response. This allows both within-pathway comparisons (PARP1 vs PARP2) and cross-pathway comparisons (BER vs Checkpoint vs NHEJ).


## Repair pathways

### BER — Base Excision Repair
PARP1 and PARP2 detect single-strand DNA breaks and recruit repair machinery. PARP inhibitors (olaparib, niraparib, rucaparib, talazoparib) are approved for BRCA-mutated breast and ovarian cancer, making this the most clinically advanced group.

### Checkpoint signaling
ATR and ATM are kinases that detect DNA damage and activate cell cycle checkpoints (pausing the cell so it can repair before dividing). ATR responds primarily to replication stress, ATM to double-strand breaks. Inhibitors of both are in clinical trials, often in combination with PARP inhibitors.

### NHEJ — Non-Homologous End Joining
DNA-PKcs is the catalytic subunit of the DNA-PK complex, which repairs double-strand breaks by directly joining the broken ends. Inhibitors are being explored in combination with radiotherapy.


## Notes

- Record counts were obtained from ChEMBL 36 using: `confidence_score >= 7`, `standard_type IN ('IC50', 'Ki')`, `standard_units = 'nM'`.
- Target IDs were verified against the `target_dictionary` table. Previous incorrect IDs (CHEMBL5882, CHEMBL4722, CHEMBL2095236) were corrected after database exploration.
