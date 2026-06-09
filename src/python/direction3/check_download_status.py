#!/usr/bin/env python3
"""Check ARCHS4 H5 download status."""
import os, sys
from pathlib import Path

H5 = Path("/home/scroll/zzhang/transcriptome-aging/data/archs4/mouse_gene_v2.5.h5")
EXPECTED = 38960132574

if not H5.exists():
    print("H5 file not found.")
    sys.exit(1)

actual = H5.stat().st_size
pct = actual / EXPECTED * 100
remaining = EXPECTED - actual
print(f"Downloaded: {actual:,} / {EXPECTED:,} bytes ({pct:.2f}%)")
print(f"Remaining:  {remaining:,} bytes ({remaining/1024**3:.1f} GB)")

# Estimate speed from log
log = Path("/home/scroll/zzhang/transcriptome-aging/logs/direction3/wget_download.log")
if log.exists():
    lines = log.read_text().strip().split("\n")
    if lines:
        last = lines[-1]
        print(f"Last log: {last.strip()}")

if actual >= EXPECTED:
    print("\nDOWNLOAD COMPLETE")
    sys.exit(0)
else:
    print("\nDownload in progress...")
    sys.exit(1)
