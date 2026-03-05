import sqlite3
import pandas as pd
from rdkit import Chem
from config import DB_PATH, TARGETS, MIN_CONFIDENCE, ACTIVITY_TYPES, ACTIVITY_UNITS

conn = sqlite3.connect(DB_PATH)

query = """
SELECT 
    md.chembl_id AS molecule_chembl_id,
    cs.canonical_smiles,
    td.chembl_id AS target_chembl_id,
    act.standard_type,
    act.standard_value,
    act.standard_units,
    act.pchembl_value
FROM activities act
JOIN assays a        ON act.assay_id = a.assay_id
JOIN target_dictionary td ON a.tid = td.tid
JOIN molecule_dictionary md ON act.molregno = md.molregno
JOIN compound_structures cs ON act.molregno = cs.molregno
WHERE a.confidence_score >= ?
AND act.standard_type IN (?, ?)
AND act.standard_units = ?
AND td.chembl_id IN (?, ?, ?, ?, ?)
"""

target_ids = list(TARGETS.keys())
params = [MIN_CONFIDENCE, *ACTIVITY_TYPES, ACTIVITY_UNITS, *target_ids]

df = pd.read_sql_query(query, conn, params=params)
print(f"Initial count: {df.count()}")

df.dropna(subset=["standard_value"], inplace=True) # remove rows where standard_value is missing
print(f"After removing rows with missing standard_value: {df.count()}")

invalid_smiles_count = 0
for index, row in df.iterrows():
    smiles = row.canonical_smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        df.drop(index, inplace=True)
        invalid_smiles_count += 1

if invalid_smiles_count > 0:
    print(f"Removed {invalid_smiles_count} rows with invalid SMILES")

# TODO: 4. STANDARDIZE molecules
#    - Strip salts (keep largest fragment)
#    - Neutralize charges
#    - Canonical tautomer

# TODO: 5. DEDUPLICATE
#    - Group by (target_chembl_id, canonical_smiles)
#    - If same molecule appears multiple times for same target:
#      keep the row with the lowest (best) IC50

# TODO: 6. COMPUTE pIC50
#    - If pchembl_value is already present -> use it
#    - If missing -> calculate: pIC50 = -log10(standard_value * 1e-9)

# TODO: 7. SAVE
#    - Save cleaned DataFrame to data/processed/cleaned_activities.csv
#    - Print summary: how many rows per target, how many removed