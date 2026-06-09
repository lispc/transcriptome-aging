#!/usr/bin/env python3
"""
GEO E-utilities search for mouse RNA-seq studies matching target queries.
Used as backup if ARCHS4 H5 download fails or for validation.
"""
import os, sys, time, json, urllib.request, urllib.parse
from pathlib import Path
from collections import defaultdict
import pandas as pd

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
OUT_DIR = PROJECT_ROOT / "data" / "geo_expansion"
OUT_DIR.mkdir(parents=True, exist_ok=True)

QUERIES = [
    ("rapamycin liver mouse RNA-seq", "rapamycin_liver"),
    ("caloric restriction liver mouse RNA-seq", "cr_liver"),
    ("metformin liver mouse RNA-seq", "metformin_liver"),
    ("aging liver mouse RNA-seq", "aging_liver"),
    ("dexamethasone mouse RNA-seq", "dexamethasone"),
    ("doxorubicin heart mouse RNA-seq", "doxorubicin_heart"),
    ("young old liver mouse RNA-seq", "young_old_liver"),
]

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def esearch(query, db="gds", retmax=500):
    url = f"{ESEARCH}?db={db}&term={urllib.parse.quote(query)}&retmax={retmax}&retmode=json"
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = json.loads(resp.read())
    return data.get("esearchresult", {}).get("idlist", [])

def esummary_ids(ids, db="gds", max_retries=3, batch_size=200):
    if not ids:
        return {}
    all_results = {}
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        idstr = ",".join(batch)
        url = f"{ESUMMARY}?db={db}&id={idstr}&retmode=json"
        for attempt in range(max_retries):
            try:
                with urllib.request.urlopen(url, timeout=120) as resp:
                    data = json.loads(resp.read())
                batch_result = data.get("result", {})
                # Merge, skipping 'uids' key
                for k, v in batch_result.items():
                    if k != "uids":
                        all_results[k] = v
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    wait = 2 ** attempt + 1
                    print(f"    Rate limited, waiting {wait}s...")
                    time.sleep(wait)
                else:
                    raise
        else:
            print(f"    Failed to fetch batch {i//batch_size + 1}")
        time.sleep(0.5)
    return all_results

def fetch_gse_summary(query_label, query_text):
    print(f"\n[{query_label}] Searching: {query_text}")
    ids = esearch(query_text)
    print(f"  Found {len(ids)} IDs")
    if not ids:
        return []
    summaries = esummary_ids(ids)
    records = []
    for uid in ids:
        info = summaries.get(str(uid), {})
        if not info:
            continue
        acc = info.get("accession", "")
        title = info.get("title", "")
        summary = info.get("summary", "")
        gds = info.get("gds", "")
        n_samples = info.get("n_samples", "")
        records.append({
            "query_label": query_label,
            "query": query_text,
            "uid": uid,
            "accession": acc,
            "title": title,
            "summary": summary,
            "gds": gds,
            "n_samples": n_samples,
        })
    print(f"  Summaries: {len(records)}")
    return records

def main():
    all_records = []
    for query_text, query_label in QUERIES:
        recs = fetch_gse_summary(query_label, query_text)
        all_records.extend(recs)
        # Save incremental results
        if recs:
            df_partial = pd.DataFrame(all_records)
            out = OUT_DIR / "geo_eutils_search_results.csv"
            df_partial.to_csv(out, index=False)
            print(f"  Incremental save: {len(df_partial)} records")
        time.sleep(1.5)  # NCBI rate limit

    df = pd.DataFrame(all_records)
    out = OUT_DIR / "geo_eutils_search_results.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved {len(df)} records to {out}")

    # Summary by query
    print("\n--- Summary ---")
    print(df.groupby("query_label").size().to_string())

if __name__ == "__main__":
    main()
