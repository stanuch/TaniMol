import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm


def fetch_activity_data(
    db_path, targets, min_confidence, activity_types, activity_units
):
    """Query ChEMBL SQLite database for bioactivity records.

    Builds a parameterized SQL query from targets dict keys and
    activity_types list. Filters by confidence_score >= min_confidence.
    Connection is opened and closed automatically.
    """
    target_ids = list(targets.keys())
    activity_placeholders = ", ".join("?" * len(activity_types))
    target_placeholders = ", ".join("?" * len(target_ids))

    query = f"""
    SELECT
        md.chembl_id  AS molecule_chembl_id,
        cs.canonical_smiles,
        td.chembl_id  AS target_chembl_id,
        act.standard_type,
        act.standard_relation,
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
      AND act.standard_relation IN ('=', '~')
      AND act.standard_units  = ?
      AND td.chembl_id       IN ({target_placeholders})
      AND act.standard_value > 0
      AND act.standard_value <= 100000
    """

    params = [min_confidence, *activity_types, activity_units, *target_ids]

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=params)

    return df


def drop_missing_values(df, column="standard_value"):
    df = df.dropna(subset=[column])
    return df


def validate_smiles(df, smiles_column="canonical_smiles"):
    """Drop rows whose SMILES cannot be parsed by RDKit (MolFromSmiles returns None)."""
    valid_mask = df[smiles_column].apply(
        lambda smi: Chem.MolFromSmiles(smi) is not None
    )
    return df[valid_mask].copy()


def standardize_molecules(df, smiles_column="canonical_smiles"):
    """Standardize SMILES in place: strip salts → neutralize → canonicalize tautomers.

    Uses RDKit's LargestFragmentChooser, Uncharger, and TautomerEnumerator
    via a performant pandas .apply() operation.
    Rows where MolFromSmiles returns None are silently skipped (SMILES unchanged).
    """
    chooser = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    tautomer_canonicalizer = rdMolStandardize.TautomerEnumerator()

    def _standardize(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return smi
        mol = chooser.choose(mol)
        mol = uncharger.uncharge(mol)
        mol = tautomer_canonicalizer.Canonicalize(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)

    df = df.copy()
    tqdm.pandas(desc="Standardizing")
    df[smiles_column] = df[smiles_column].progress_apply(_standardize)
    return df


def deduplicate(
    df,
    target_column="target_chembl_id",
    smiles_column="canonical_smiles",
    value_column="standard_value",
):
    """Keep one row per (target, SMILES) pair taking the geometric median IC50.

    Since IC50 is log-normally distributed, the geometric median is more
    statistically sound than the arithmetic median. Other columns keep
    their first occurrence.
    """
    # Znajdź tylko poprawne (dodatnie) wartości, żeby móc zlogarytmować
    valid_mask = df[value_column] > 0

    medians = (
        df[valid_mask]
        .groupby([target_column, smiles_column])[value_column]
        .apply(lambda x: np.exp(np.median(np.log(x))))
        .reset_index()
    )

    df_dedup = df.drop_duplicates(
        subset=[target_column, smiles_column], keep="first"
    ).copy()
    df_dedup = df_dedup.drop(columns=[value_column])

    df = pd.merge(df_dedup, medians, on=[target_column, smiles_column], how="left")

    # Remove compounds that had ONLY non-positive values (they became NaN after merge)
    df = df.dropna(subset=[value_column])

    return df


def compute_pic50(df, pchembl_column="pchembl_value", value_column="standard_value"):
    """Fill missing pchembl_value using pIC50 = -log10(IC50_nM x 1e-9).

    Existing pchembl values are kept as-is. Rows with IC50 ≤ 0 will
    remain NaN (log of non-positive is undefined) — drop them afterwards.
    Expects standard_value in nM.
    """
    for i, row in tqdm(df.iterrows(), total=len(df), desc="Computing pIC50"):
        if pd.isna(df.loc[i, pchembl_column]):
            sv = df.loc[i, value_column]
            if sv is not None and sv > 0:
                df.loc[i, pchembl_column] = -np.log10(sv * 1e-9)
    return df


def save_cleaned_data(df, output_path, targets=None):
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} rows to {output_path}")
    if targets:
        print("\nRows per target:")
        for chembl_id, info in targets.items():
            count = (df["target_chembl_id"] == chembl_id).sum()
            print(f"- {info['name']} ({chembl_id}): {count}")
