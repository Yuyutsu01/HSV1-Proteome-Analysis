#!/usr/bin/env python3
"""
Setup Script: Populate HSV_Proteome_Analysis_Workspace at root
"""

import os
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
WORKSPACE = ROOT / "HSV_Proteome_Analysis_Workspace"
SRC = ROOT / "HSV_Computational_analysis"

def main():
    print(f"Populating Workspace at: {WORKSPACE}")
    
    # 1. Raw Data
    raw_dst = WORKSPACE / "data" / "raw"
    raw_dst.mkdir(parents=True, exist_ok=True)
    for f in (SRC / "data" / "raw").glob("*"):
        if f.is_file():
            shutil.copy2(f, raw_dst / f.name)
            
    # 2. Intermediate Data
    inter_dst = WORKSPACE / "data" / "intermediate"
    inter_dst.mkdir(parents=True, exist_ok=True)
    for f in (SRC / "data" / "intermediate").glob("*"):
        if f.is_file():
            shutil.copy2(f, inter_dst / f.name)
            
    # 3. Processed Data
    proc_dst = WORKSPACE / "data" / "processed"
    proc_dst.mkdir(parents=True, exist_ok=True)
    for f in (SRC / "data" / "processed").glob("*"):
        if f.is_file():
            shutil.copy2(f, proc_dst / f.name)
            
    # Copy embeddings
    esm2_dst = proc_dst / "embeddings" / "esm2"
    esm2_dst.mkdir(parents=True, exist_ok=True)
    for f in (SRC / "data" / "processed" / "embeddings" / "esm2").glob("*"):
        if f.is_file():
            shutil.copy2(f, esm2_dst / f.name)
            
    # Copy kmers
    kmer_dst = proc_dst / "kmer_features"
    kmer_dst.mkdir(parents=True, exist_ok=True)
    for f in (SRC / "data" / "processed" / "phase3_kmer_features").glob("*"):
        if f.is_file():
            shutil.copy2(f, kmer_dst / f.name)

    # 4. Tables
    for sub, prefix in [("phase1", "phase1"), ("phase2", "phase2"), ("phase3", "phase3"), ("phase4", "phase4"), ("phase4_1", "phase4_1")]:
        dst = WORKSPACE / "tables" / sub
        dst.mkdir(parents=True, exist_ok=True)
        for f in (SRC / "results" / "tables").glob(f"{prefix}*.csv"):
            shutil.copy2(f, dst / f.name)
            
    # Also copy final_ tables into tables/phase2
    for f in (SRC / "results" / "tables").glob("final_*.csv"):
        shutil.copy2(f, WORKSPACE / "tables" / "phase2" / f.name)

    # 5. Figures
    fig3_dst = WORKSPACE / "figures" / "phase3"
    fig3_dst.mkdir(parents=True, exist_ok=True)
    for f in (SRC / "results" / "figures" / "phase3").glob("*.png"):
        shutil.copy2(f, fig3_dst / f.name)

    fig4_dst = WORKSPACE / "figures" / "phase4"
    fig4_dst.mkdir(parents=True, exist_ok=True)
    for f in (SRC / "results" / "figures" / "phase4").glob("*.png"):
        shutil.copy2(f, fig4_dst / f.name)

    print("All files successfully copied into workspace.")

if __name__ == "__main__":
    main()
