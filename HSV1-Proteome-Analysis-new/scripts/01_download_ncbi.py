"""
Script 01: Programmatic NCBI GenBank Record Acquisition

Scientific Objective:
Download the official reference genome and proteome record for Herpes Simplex Virus Type 1
(Human herpesvirus 1 strain 17, accession: NC_001806.2) programmatically from NCBI Entrez.
Maintains full cryptographic and provenance tracking (SHA256 checksum, timestamp, URL, versioning).

Author: Computational Biology Pipeline
"""

import os
import sys
import json
import hashlib
import datetime
import urllib.request
import yaml
from Bio import Entrez, __version__ as biopython_version

def load_config(config_path="config.yaml"):
    """Load pipeline YAML configuration."""
    if not os.path.exists(config_path):
        # Check parent or current directory
        if os.path.exists(os.path.join("..", config_path)):
            config_path = os.path.join("..", config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def compute_sha256(filepath):
    """Compute SHA256 checksum of a file for cryptographic provenance."""
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            hasher.update(chunk)
    return hasher.hexdigest()

def download_ncbi_genbank(config_path="config.yaml"):
    """
    Download reference GenBank flatfile from NCBI Entrez with complete provenance tracking.
    """
    cfg = load_config(config_path)
    accession = cfg['ncbi']['accession']
    db = cfg['ncbi']['database']
    rettype = cfg['ncbi']['rettype']
    retmode = cfg['ncbi']['retmode']
    email = os.environ.get("NCBI_EMAIL", cfg['ncbi'].get('email', 'researcher@example.com'))
    
    out_dir = cfg['paths']['data_raw_ncbi']
    os.makedirs(out_dir, exist_ok=True)
    target_path = os.path.join(out_dir, f"{accession}.gbk")
    manifest_path = os.path.join(out_dir, "download_manifest.json")
    
    print(f"[Phase 1] Initializing NCBI Entrez fetch for Accession: {accession}")
    print(f"  - Database: {db}")
    print(f"  - Email configured: {email}")
    print(f"  - Target file: {target_path}")
    
    # Configure Entrez
    Entrez.email = email
    download_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={db}&id={accession}&rettype={rettype}&retmode={retmode}"
    
    # Perform download if not already cached or verify existing
    if os.path.exists(target_path) and os.path.getsize(target_path) > 10000:
        print(f"[Phase 1] Local file already exists ({os.path.getsize(target_path)} bytes). Verifying...")
    else:
        print(f"[Phase 1] Downloading from NCBI E-Utilities: {download_url}")
        req = urllib.request.Request(download_url, headers={'User-Agent': f'BioPython/{biopython_version} (NCBI_Downloader)'})
        with urllib.request.urlopen(req) as resp, open(target_path, "wb") as out_f:
            out_f.write(resp.read())
        print(f"[Phase 1] Download complete: {target_path} ({os.path.getsize(target_path)} bytes)")
        
    checksum = compute_sha256(target_path)
    manifest = {
        "accession": accession,
        "database": db,
        "rettype": rettype,
        "retmode": retmode,
        "download_timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "source_url": download_url,
        "local_file_path": target_path,
        "file_size_bytes": os.path.getsize(target_path),
        "sha256_checksum": checksum,
        "biopython_version": biopython_version,
        "python_version": sys.version
    }
    
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
        
    print(f"[Phase 1] Provenance manifest saved to {manifest_path}")
    print(f"  - SHA256: {checksum}")
    return target_path, manifest

if __name__ == "__main__":
    download_ncbi_genbank()
