"""
Geneformer Tokenization for Tabula Muris Senis data.
Convert processed scRNA-seq data to rank-value encoding for Geneformer.
"""
import os
import sys
from pathlib import Path
import scanpy as sc

# Add Geneformer to path
sys.path.insert(0, "/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2")

from geneformer import TranscriptomeTokenizer

DATA_DIR = Path("/home/scroll/zzhang/transcriptome-aging/data/scrna_aging/tabula_muris_senis")
OUTPUT_DIR = Path("/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def tokenize_liver():
    """Tokenize TMS Liver data for Geneformer."""
    input_path = DATA_DIR / "tms_liver_10x.h5ad"
    
    print(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    print(f"Loaded: {adata.shape}")
    
    # Geneformer expects specific column names
    # Ensure required metadata columns exist
    if 'age_group' not in adata.obs.columns:
        # Recompute age_group
        def categorize_age(age_str):
            try:
                age_num = float(str(age_str).replace('m', '').replace('M', ''))
                if age_num <= 6:
                    return 'young'
                elif age_num <= 12:
                    return 'middle'
                else:
                    return 'old'
            except:
                return 'unknown'
        adata.obs['age_group'] = adata.obs['age'].apply(categorize_age)
    
    # For Geneformer, we need Ensembl gene IDs.
    # TMS uses gene symbols. We need to map them.
    # For now, use symbols as-is (Geneformer vocab includes many symbols)
    print(f"Gene names (first 10): {list(adata.var_names[:10])}")
    
    # Check overlap with Geneformer vocab
    import json
    vocab_path = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2/geneformer/token_dictionary.json"
    if os.path.exists(vocab_path):
        with open(vocab_path) as f:
            vocab = json.load(f)
        overlap = sum(1 for g in adata.var_names if g in vocab)
        print(f"Genes in Geneformer vocab: {overlap} / {len(adata.var_names)} ({100*overlap/len(adata.var_names):.1f}%)")
    
    # Save with required metadata
    processed_path = OUTPUT_DIR / "tms_liver_for_tokenization.h5ad"
    adata.write_h5ad(processed_path)
    print(f"Saved pre-tokenization data to {processed_path}")
    
    # Tokenize
    print("\nStarting tokenization...")
    tk = TranscriptomeTokenizer(
        nproc=4,
        model_input_size=4096,  # V2-104M uses 4096
        # Geneformer V2 expects specific column for tokenization
    )
    
    # Tokenize with age_group as the key metadata
    tk.tokenize_data(
        adata,
        output_directory=str(OUTPUT_DIR),
        output_prefix="tms_liver_tokenized",
        file_format="h5ad",
        use_generator=True,
    )
    
    print(f"\nTokenization complete! Output in {OUTPUT_DIR}")


if __name__ == "__main__":
    tokenize_liver()
