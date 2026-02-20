# === ChEMBL database version ===

CHEMBL_VERSION = "36"
DB_FILENAME = f"chembl_{CHEMBL_VERSION}.db"


# === Targets for the database queries ===

TARGETS = {
    "CHEMBL3105": {"name": "PARP1", "pathway": "BER"},
    "CHEMBL5366": {"name": "PARP2", "pathway": "BER"},
    "CHEMBL5024": {"name": "ATR",   "pathway": "Checkpoint"},
    "CHEMBL3797": {"name": "ATM",   "pathway": "Checkpoint"},
    "CHEMBL3142": {"name": "DNA-PKcs", "pathway": "NHEJ"},
}


# === Filters ===

MIN_CONFIDENCE = 7
ACTIVITY_TYPES = ("IC50", "Ki")
ACTIVITY_UNITS = "nM"