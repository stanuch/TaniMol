# ===========================
# ======= PSUEDO CODE =======
# ===========================

# INPUT:  path to chembl_36.db
# OUTPUT: cleaned DataFrame saved as CSV to data/processed/

# 1. CONNECT to SQLite database

# 2. QUERY activities for all targets in TARGETS dict
#    - JOIN activities + assays + target_dictionary + molecule_dictionary + compound_structures
#    - WHERE confidence >= 7, type IN (IC50, Ki), units = nM
#    - Store result as DataFrame

# 3. VALIDATE SMILES
#    - For each row, try to parse canonical_smiles with RDKit
#    - Remove rows where SMILES is invalid (RDKit returns None)
#    - Log how many were removed

# 4. STANDARDIZE molecules
#    - Strip salts (keep largest fragment)
#    - Neutralize charges
#    - Canonical tautomer

# 5. DEDUPLICATE
#    - Group by (target_chembl_id, canonical_smiles)
#    - If same molecule appears multiple times for same target:
#      keep the row with the lowest (best) IC50

# 6. COMPUTE pIC50
#    - If pchembl_value is already present -> use it
#    - If missing -> calculate: pIC50 = -log10(standard_value * 1e-9)

# 7. SAVE
#    - Save cleaned DataFrame to data/processed/cleaned_activities.csv
#    - Print summary: how many rows per target, how many removed