"""Fingerprint generation (Morgan/ECFP, MACCS, RDKit topological) for Tanimoto similarity."""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, RDKFingerprint
from tqdm import tqdm


def mol_from_smiles(smiles):
    """Parse SMILES to an RDKit Mol. Returns None if SMILES is invalid or empty."""
    raise NotImplementedError


def generate_morgan_fp(mol, radius=2, n_bits=2048):
    """Generate Morgan (ECFP) fingerprint as a NumPy bit array.

    radius=2 corresponds to ECFP4, radius=3 to ECFP6.
    mol must be a sanitized RDKit Mol (not None).

    - Use AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    - Convert to NumPy: np.array(fp)
    """
    raise NotImplementedError


def generate_maccs_fp(mol):
    """Generate MACCS keys (166 predefined substructure patterns).

    Returns an array of shape (167,) — bit 0 is always unused.
    mol must be a sanitized RDKit Mol (not None).

    - Use MACCSkeys.GenMACCSKeys(mol)
    - Convert to NumPy: np.array(fp)
    """
    raise NotImplementedError


def generate_rdkit_fp(mol, fp_size=2048):
    """Generate RDKit topological fingerprint (path-based, not circular).
    mol must be a sanitized RDKit Mol (not None).

    - Use Chem.RDKFingerprint(mol, fpSize=fp_size)
    - Convert to NumPy: np.array(fp)
    """
    raise NotImplementedError


def add_fingerprints(df, smiles_column="canonical_smiles", fp_type="morgan", **kwargs):
    """Generate fingerprints for all molecules in a DataFrame.

    fp_type must be 'morgan', 'maccs', or 'rdkit'. Extra **kwargs are
    forwarded to the generator (e.g. radius, n_bits for Morgan).
    Molecules that fail parsing are appended as None.

    1. Pick the right generator function based on fp_type
    2. Loop through df[smiles_column] with tqdm
    3. For each SMILES: mol_from_smiles → generate fp
    4. Handle None mols gracefully (append None or skip)
    5. Return the list of fingerprints
    """
    raise NotImplementedError
