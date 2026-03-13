"""
Preprocessing pipeline for ChEMBL activity data.

Functions for fetching, cleaning, and preparing bioactivity data
for downstream analysis in notebooks.
"""

import sqlite3
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm


def fetch_activity_data(db_path, targets, min_confidence, activity_types, activity_units):
    """
    Fetch activity data from ChEMBL SQLite database.

    Parameters
    ----------
    db_path : str or Path
        Path to the ChEMBL SQLite database file.
    targets : dict
        Dictionary with target ChEMBL IDs as keys.
    min_confidence : int
        Minimum assay confidence score.
    activity_types : list of str
        Activity types to include (e.g. ["IC50", "Ki"]).
    activity_units : str
        Expected unit for standard_value (e.g. "nM").

    Returns
    -------
    pd.DataFrame
        Raw activity data.
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

    params = [min_confidence, *activity_types, activity_units, *target_ids]

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(query, conn, params=params)

    return df


def drop_missing_values(df, column="standard_value"):
    """
    Drop rows where the given column has missing values.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    column : str
        Column to check for missing values.

    Returns
    -------
    pd.DataFrame
        DataFrame with missing values removed.
    """
    df = df.dropna(subset=[column])
    return df


def validate_smiles(df, smiles_column="canonical_smiles"):
    """
    Validate SMILES strings using RDKit and remove invalid ones.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    smiles_column : str
        Name of the column containing SMILES strings.

    Returns
    -------
    pd.DataFrame
        DataFrame with only valid SMILES.
    """
    valid_mask = df[smiles_column].apply(
        lambda smi: Chem.MolFromSmiles(smi) is not None
    )
    return df[valid_mask].copy()


def standardize_molecules(df, smiles_column="canonical_smiles"):
    """
    Standardize molecular structures.

    Steps:
    - Strip salts (keep largest fragment)
    - Neutralize charges
    - Canonical tautomer

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    smiles_column : str
        Name of the column containing SMILES strings.

    Returns
    -------
    pd.DataFrame
        DataFrame with standardized SMILES.
    """
    chooser = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    tautomer_canonicalizer = rdMolStandardize.TautomerEnumerator()

    for i, row in tqdm(df.iterrows(), total=len(df), desc="Standardizing"):
        mol = Chem.MolFromSmiles(row[smiles_column])
        if mol is None:
            continue
        mol = chooser.choose(mol)
        mol = uncharger.uncharge(mol)
        mol = tautomer_canonicalizer.Canonicalize(mol)
        df.loc[i, smiles_column] = Chem.MolToSmiles(mol)
    return df


def deduplicate(df, target_column="target_chembl_id", smiles_column="canonical_smiles",
                value_column="standard_value"):
    """
    Deduplicate by (target, SMILES) keeping the row with the lowest (best) IC50.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    target_column : str
        Column with target identifiers.
    smiles_column : str
        Column with SMILES strings.
    value_column : str
        Column with activity value (lower = better).

    Returns
    -------
    pd.DataFrame
        Deduplicated DataFrame.
    """
    df = df.sort_values(by=value_column)
    df = df.drop_duplicates(subset=[target_column, smiles_column], keep="first")
    return df


def compute_pic50(df, pchembl_column="pchembl_value", value_column="standard_value"):
    """
    Compute pIC50 values.

    If pchembl_value is already present, use it.
    If missing, calculate: pIC50 = -log10(standard_value * 1e-9).

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    pchembl_column : str
        Column with existing pChEMBL values.
    value_column : str
        Column with standard_value in nM.

    Returns
    -------
    pd.DataFrame
        DataFrame with a 'pic50' column.
    """
    for i, row in tqdm(df.iterrows(), total=len(df), desc="Computing pIC50"):
        if pd.isna(df.loc[i, pchembl_column]):
            sv = df.loc[i, value_column]
            if sv is not None and sv > 0:
                df.loc[i, pchembl_column] = -np.log10(sv * 1e-9)
    return df


def save_cleaned_data(df, output_path, targets=None):
    """
    Save cleaned DataFrame to CSV and print summary.

    Parameters
    ----------
    df : pd.DataFrame
        Cleaned DataFrame.
    output_path : str or Path
        Path for the output CSV file.
    targets : dict, optional
        Target dictionary for summary reporting.
    """
    df.to_csv(output_path, index=False)
    print(f"Saved {len(df)} rows to {output_path}")
    if targets:
        print("\nRows per target:")
        for chembl_id, info in targets.items():
            count = (df["target_chembl_id"] == chembl_id).sum()
            print(f"- {info['name']} ({chembl_id}): {count}")
