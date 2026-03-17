from pathlib import Path
from datetime import datetime

# === ChEMBL database version ===

CHEMBL_VERSION = "36"
DB_FILENAME = f"chembl_{CHEMBL_VERSION}.db"


# === Output filename (change if needed) ===

OUTPUT_FILE_NAME = "cleaned_activities.csv"


# === Filters ===

MIN_CONFIDENCE = 7
ACTIVITY_TYPES = ["IC50"]
ACTIVITY_UNITS = "nM"
MORGAN_RADIUS = 2
CLUSTERING_THRESHOLD = 0.6


# === Targets for the database queries ===

TARGETS = {
    "CHEMBL3105": {"name": "PARP1", "pathway": "BER"},
    "CHEMBL5366": {"name": "PARP2", "pathway": "BER"},
    "CHEMBL5024": {"name": "ATR",   "pathway": "Checkpoint"},
    "CHEMBL3797": {"name": "ATM",   "pathway": "Checkpoint"},
    "CHEMBL3142": {"name": "DNA-PKcs", "pathway": "NHEJ"},
}


# === Paths ===

# Project Root
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
OUTPUT_PATH = PROCESSED_DIR / OUTPUT_FILE_NAME

# Fingerprints
FINGERPRINTS_DIR = PROCESSED_DIR / "fingerprints"
MORGAN_FP_PATH = FINGERPRINTS_DIR / "morgan_fps.npy"
MACCS_FP_PATH = FINGERPRINTS_DIR / "maccs_fps.npy"
RDKIT_FP_PATH = FINGERPRINTS_DIR / "rdkit_fps.npy"

# Similarity Matrix
SIMILARITY_MATRIX_DIR = PROCESSED_DIR / "similarity"
MORGAN_SIM_PATH = SIMILARITY_MATRIX_DIR / "morgan_sim.npy"
MACCS_SIM_PATH = SIMILARITY_MATRIX_DIR / "maccs_sim.npy"
RDKIT_SIM_PATH = SIMILARITY_MATRIX_DIR / "rdkit_sim.npy"

# Database
URL = f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_{CHEMBL_VERSION}_sqlite.tar.gz"
ARCHIVE_PATH = SRC_DIR / f"chembl_{CHEMBL_VERSION}_sqlite.tar.gz"
DB_PATH = RAW_DIR / f"chembl_{CHEMBL_VERSION}.db"
