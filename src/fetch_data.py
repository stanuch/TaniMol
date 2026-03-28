import os
import sys
import shutil
import glob
import urllib.request
import tarfile

# Configuration (version, targets, filters) is in config.py
# To change the ChEMBL database version, edit CHEMBL_VERSION in config.py

from config import CHEMBL_VERSION, URL, ARCHIVE_PATH, DB_PATH, RAW_DIR, SRC_DIR

if os.path.exists(DB_PATH):
    print(f"Database already exists at {DB_PATH}, skipping downloading process.")
else:
    if os.path.exists(ARCHIVE_PATH):
        print("Archive already downloaded, skipping download.")
    else:
        print(f"ChEMBL {CHEMBL_VERSION} database not found at {DB_PATH}")

        while True:
            answer = input("Do you want to start downloading? [Y/n]: ").strip().lower()
            if answer in ("y", "n"):
                break
            print("Please enter 'y' or 'n'.")
        if answer == "n":
            print("Download cancelled.")
            sys.exit(0)
        print(f"\nDownloading ChEMBL Database, Version: {CHEMBL_VERSION}...")

        def progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            if total_size > 0:
                percent = min(downloaded / total_size * 100, 100)
                mb_down = downloaded / (1024 * 1024)
                mb_total = total_size / (1024 * 1024)
                filled = int(40 * percent / 100)
                bar = "█" * filled + "_" * (40 - filled)
                sys.stdout.write(
                    f"\r  [{bar}] {percent:5.1f}%  ({mb_down:.0f}/{mb_total:.0f} MB)"
                )
                sys.stdout.flush()

        urllib.request.urlretrieve(URL, ARCHIVE_PATH, reporthook=progress)
        print("\nDownloaded!")

    print("Extracting... (May take a few minutes)")
    with tarfile.open(ARCHIVE_PATH, "r:gz") as tar:
        tar.extractall(path=SRC_DIR)
    print("Extracted!")

    db_files = glob.glob(os.path.join(SRC_DIR, "**", "*.db"), recursive=True)

    if not db_files:
        print("ERROR: Could not find .db file after extraction!")
    else:
        os.makedirs(RAW_DIR, exist_ok=True)
        shutil.move(db_files[0], DB_PATH)
        print(f"Moved database to {DB_PATH}")

    if os.path.exists(ARCHIVE_PATH):
        os.remove(ARCHIVE_PATH)
        print("Removed .tar.gz archive.")

    extracted_folder = os.path.join(SRC_DIR, f"chembl_{CHEMBL_VERSION}")
    if os.path.exists(extracted_folder):
        shutil.rmtree(extracted_folder)
        print("Removed leftover extraction folders.")

    print(f"Done! Database is at: {DB_PATH}")
