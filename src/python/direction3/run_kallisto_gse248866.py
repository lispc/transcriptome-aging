#!/usr/bin/env python3
import subprocess
from pathlib import Path

PROJECT = Path(__file__).parent.parent.parent.parent
IDX = PROJECT / "data" / "ref" / "mus_musculus" / "transcriptome_v52.idx"
FQDIR = PROJECT / "data" / "gse248866_fastq"
OUTDIR = PROJECT / "data" / "gse248866_kallisto"

samples = {
    "AL_1": ("SRR26975750",),
    "AL_2": ("SRR26975749",),
    "AL_3": ("SRR26975748",),
    "CR_1": ("SRR26975747",),
    "CR_2": ("SRR26975746",),
    "CR_3": ("SRR26975745",),
    "NR_1": ("SRR33542637",),
    "NR_2": ("SRR33542636",),
    "NR_3": ("SRR33542635",),
    "Cort_1": ("SRR33542634",),
    "Cort_2": ("SRR33542633",),
    "Cort_3": ("SRR33542632",),
}

for name, (srr,) in samples.items():
    out = OUTDIR / name
    if (out / "abundance.tsv").exists():
        print(f"SKIP {name}")
        continue
    
    r1 = FQDIR / f"{srr}_1.fastq.gz"
    r2 = FQDIR / f"{srr}_2.fastq.gz"
    if not r1.exists():
        r1 = FQDIR / f"{srr}_1.fastq"
    if not r2.exists():
        r2 = FQDIR / f"{srr}_2.fastq"
    
    print(f"QUANT {name}: {r1.name} + {r2.name}")
    subprocess.run([
        "kallisto", "quant", "-i", str(IDX), "-o", str(out),
        "--threads", "4", str(r1), str(r2)
    ], check=True)
    print(f"  DONE {name}")

print("\nAll quantification complete!")
