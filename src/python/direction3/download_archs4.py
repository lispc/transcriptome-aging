#!/usr/bin/env python3
"""
Download ARCHS4 mouse_gene_v2.5.h5 with resume support.
"""
import os, sys, urllib.request, time

URL = "https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.5.h5"
OUT = "data/archs4/mouse_gene_v2.5.h5"
BLOCK = 1024 * 1024  # 1 MB

def download():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    existing = os.path.getsize(OUT) if os.path.exists(OUT) else 0
    req = urllib.request.Request(URL)
    if existing:
        req.add_header("Range", f"bytes={existing}-")
        print(f"Resuming from {existing} bytes")
    else:
        print("Starting fresh download")

    t0 = time.time()
    with urllib.request.urlopen(req) as resp:
        total = resp.headers.get("Content-Length")
        total = int(total) + existing if total else None
        mode = "ab" if existing else "wb"
        with open(OUT, mode) as f:
            downloaded = existing
            while True:
                chunk = resp.read(BLOCK)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if time.time() - t0 > 5:
                    pct = f"{downloaded/total*100:.1f}%" if total else f"{downloaded}"
                    print(f"  {pct} ({downloaded}/{total}) bytes")
                    t0 = time.time()
    print(f"Done: {OUT} ({os.path.getsize(OUT)} bytes)")

if __name__ == "__main__":
    download()
