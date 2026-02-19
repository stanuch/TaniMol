import os
import sys
import shutil
import glob
import urllib.request
import tarfile

# Paths

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")

# Configuration

CHEMBL_VERSION = "36"
URL = f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_{CHEMBL_VERSION}_sqlite.tar.gz"

ARCHIVE_PATH = os.path.join(SCRIPT_DIR, f"chembl_{CHEMBL_VERSION}_sqlite.tar.gz")
DB_PATH = os.path.join(RAW_DIR, f"chembl_{CHEMBL_VERSION}.db")

if os.path.exists(DB_PATH):
    print(f"Database already exists at {DB_PATH}, skipping downloading process.")
else:
    if os.path.exists(ARCHIVE_PATH):
        print("Archive already downloaded, skipping.")
    else:
        print(f"Downloading ChEMBL Database, Version: {CHEMBL_VERSION}...")

        def progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            if total_size > 0:
                percent = min(downloaded / total_size * 100, 100)
                mb_down = downloaded / (1024 * 1024)
                mb_total = total_size / (1024 * 1024)
                filled = int(40 * percent / 100)
                bar = "█" * filled + "_" * (40 - filled)
                sys.stdout.write(f"\r  [{bar}] {percent:5.1f}%  ({mb_down:.0f}/{mb_total:.0f} MB)")
                sys.stdout.flush()

        urllib.request.urlretrieve(URL, ARCHIVE_PATH, reporthook=progress)
        print("\nDownloaded!")

    print("Extracting... (May take a few minutes)")
    with tarfile.open(ARCHIVE_PATH, "r:gz") as tar:
        tar.extractall(path=SCRIPT_DIR)
    print("Extracted!")

    db_files = glob.glob(os.path.join(SCRIPT_DIR, "**", "*.db"), recursive=True)

    if not db_files:
        print("ERROR: Could not find .db file after extraction!")
    else:
        os.makedirs(RAW_DIR, exist_ok=True)
        shutil.move(db_files[0], DB_PATH)
        print(f"Moved database to {DB_PATH}")

    if os.path.exists(ARCHIVE_PATH):
        os.remove(ARCHIVE_PATH)
        print("Removed .tar.gz archive.")

    extracted_folder = os.path.join(SCRIPT_DIR, f"chembl_{CHEMBL_VERSION}")
    if os.path.exists(extracted_folder):
        shutil.rmtree(extracted_folder)
        print("Removed leftover extraction folders.")

    print(f"Done! Database is at: {DB_PATH}")