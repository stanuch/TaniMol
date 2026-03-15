"""Fingerprint generation (Morgan/ECFP, MACCS, RDKit topological) for Tanimoto similarity."""

import numpy as np
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, RDKFingerprint
from tqdm import tqdm


def mol_from_smiles(smiles):
    """Parse SMILES to an RDKit Mol. Returns None if SMILES is invalid or empty."""
    if not smiles:
        return None
    return Chem.MolFromSmiles(smiles)



def generate_morgan_fp(mol, radius=2, n_bits=2048):
    """Generate Morgan (ECFP) fingerprint as a NumPy bit array.

    radius=2 corresponds to ECFP4, radius=3 to ECFP6.
    mol must be a sanitized RDKit Mol (not None).
    """
    generator = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fp = generator.GetFingerprintAsNumPy(mol)
    return fp


def generate_maccs_fp(mol):
    """Generate MACCS keys (166 predefined substructure patterns).

    Returns an array of shape (167,) — bit 0 is always unused.
    mol must be a sanitized RDKit Mol (not None).
    """
    fp = MACCSkeys.GenMACCSKeys(mol)
    return np.array(fp)


def generate_rdkit_fp(mol, fp_size=2048):
    """Generate RDKit topological fingerprint (path-based, not circular).
    mol must be a sanitized RDKit Mol (not None).
    """
    fp = Chem.RDKFingerprint(mol, fpSize=fp_size)
    return np.array(fp)


def add_fingerprints(df, fp_type, smiles_column="canonical_smiles", **kwargs):
    """Generate fingerprints for all molecules in a DataFrame.

    fp_type must be 'morgan', 'maccs', or 'rdkit'. Extra **kwargs are
    forwarded to the generator (e.g. radius, n_bits for Morgan).
    Molecules that fail parsing are appended as None.
    """
    if fp_type == "morgan":
        generator = generate_morgan_fp
    elif fp_type == "maccs":
        generator = generate_maccs_fp
    elif fp_type == "rdkit":
        generator = generate_rdkit_fp
    else:
        raise ValueError("fp_type must be 'morgan', 'maccs', or 'rdkit'")

    fingerprints = []

    for smiles in tqdm(df[smiles_column]):
        mol = mol_from_smiles(smiles)
        if mol is not None:
            fp = generator(mol, **kwargs)
            fingerprints.append(fp)
        else:
            fingerprints.append(None)

    return fingerprints


def save_fingerprints(fps, path):

    none_count = sum(1 for fp in fps if fp is None)
    if none_count > 0:
        raise ValueError(
            f"{none_count}/{len(fps)} fingerprints are None — "
            "clean your data before saving."
        )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, np.array(fps))
    print(f"Saved {len(fps)} fingerprints to {path}")


def load_fingerprints(path):
    fingerprints = np.load(path)
    return fingerprints
