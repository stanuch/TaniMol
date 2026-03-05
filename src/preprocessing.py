import sqlite3
import pandas as pd
from rdkit import Chem
from config import (
    DB_PATH,
    TARGETS,
    MIN_CONFIDENCE,
    ACTIVITY_TYPES,
    ACTIVITY_UNITS,
)

query = """
SELECT
    md.chembl_id  AS molecule_chembl_id,
    cs.canonical_smiles,
    td.chembl_id  AS target_chembl_id,
    act.standard_type,
    act.standard_value,
    act.standard_units,
    act.pchembl_value
FROM activities act
JOIN assays              a  ON act.assay_id = a.assay_id
JOIN target_dictionary   td ON a.tid        = td.tid
JOIN molecule_dictionary md ON act.molregno = md.molregno
JOIN compound_structures cs ON act.molregno = cs.molregno
WHERE a.confidence_score >= ?
  AND act.standard_type  IN ({activity_placeholders})
  AND act.standard_units  = ?
  AND td.chembl_id       IN ({target_placeholders})
"""

target_ids = list(TARGETS.keys())
activity_placeholders = ", ".join("?" * len(ACTIVITY_TYPES))
target_placeholders = ", ".join("?" * len(target_ids))
query = query.format(
    activity_placeholders=activity_placeholders,
    target_placeholders=target_placeholders,
)
params = [MIN_CONFIDENCE, *ACTIVITY_TYPES, ACTIVITY_UNITS, *target_ids]

with sqlite3.connect(DB_PATH) as conn:
    df = pd.read_sql_query(query, conn, params=params)

print(f"Initial count:\n{df.count()}")

df.dropna(subset=["standard_value"], inplace=True) # Drop rows with missing values
print(f"\nAfter removing rows with missing standard_value:\n{df.count()}")

# Validate SMILES

for index, row in df.iterrows():
    smiles = row.canonical_smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        df.drop(index, inplace=True)

# TODO: 4. Standardize molecules
#    - Strip salts (keep largest fragment)
#    - Neutralize charges
#    - Canonical tautomer

# TODO: 5. Deduplicate
#    - Group by (target_chembl_id, canonical_smiles)
#    - If same molecule appears multiple times for same target:
#      keep the row with the lowest (best) IC50

# TODO: 6. Compute pIC50
#    - If pchembl_value is already present -> use it
#    - If missing -> calculate: pIC50 = -log10(standard_value * 1e-9)

# TODO: 7. Save
#    - Save cleaned DataFrame to data/processed/cleaned_activities.csv
#    - Print summary: how many rows per target, how many removed