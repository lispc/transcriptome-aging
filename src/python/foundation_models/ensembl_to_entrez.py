"""
Map human Ensembl gene IDs to Entrez Gene IDs for tAge compatibility.
Uses mygene.info API.
"""
import json
import pickle
import mygene
from pathlib import Path

DATA_DIR = Path("/home/scroll/zzhang/transcriptome-aging/data")
OUTPUT_DIR = DATA_DIR / "gene_id_mappings"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def get_geneformer_ensembl_ids():
    """Get all Ensembl IDs from Geneformer vocab."""
    vocab_path = DATA_DIR / "geneformer/weights/geneformer-v2/geneformer/token_dictionary_gc104M.pkl"
    with open(vocab_path, 'rb') as f:
        vocab = pickle.load(f)
    
    # Filter out special tokens, keep only ENSG
    ensg_ids = [k for k in vocab.keys() if isinstance(k, str) and k.startswith('ENSG')]
    print(f"Geneformer vocab: {len(vocab)} total, {len(ensg_ids)} ENSG IDs")
    return ensg_ids


def map_ensembl_to_entrez(ensg_list, batch_size=1000):
    """Map Ensembl IDs to Entrez IDs using mygene.info."""
    mg = mygene.MyGeneInfo()
    
    results = []
    for i in range(0, len(ensg_list), batch_size):
        batch = ensg_list[i:i+batch_size]
        print(f"Mapping batch {i//batch_size + 1}/{(len(ensg_list)-1)//batch_size + 1} ({len(batch)} genes)...")
        
        query = mg.querymany(
            batch,
            scopes='ensembl.gene',
            fields='entrezgene',
            species='human',
            as_dataframe=True,
            verbose=False
        )
        results.append(query)
    
    import pandas as pd
    df = pd.concat(results)
    
    # Clean up
    df = df.reset_index()
    df = df[['query', 'entrezgene']].drop_duplicates(subset=['query'])
    df.columns = ['ensembl_id', 'entrez_id']
    df = df.dropna(subset=['entrez_id'])
    df['entrez_id'] = df['entrez_id'].astype(int).astype(str)
    
    return df


def main():
    print("=" * 60)
    print("Gene ID Mapping: Ensembl → Entrez")
    print("=" * 60)
    
    # Get Geneformer ENSG IDs
    ensg_ids = get_geneformer_ensembl_ids()
    
    # Map to Entrez
    mapping_df = map_ensembl_to_entrez(ensg_ids)
    
    # Save
    output_path = OUTPUT_DIR / "ensembl_to_entrez_geneformer.json"
    mapping_dict = dict(zip(mapping_df['ensembl_id'], mapping_df['entrez_id']))
    with open(output_path, 'w') as f:
        json.dump(mapping_dict, f)
    
    print(f"\nMapping complete!")
    print(f"  Total ENSG: {len(ensg_ids)}")
    print(f"  Mapped to Entrez: {len(mapping_dict)} ({100*len(mapping_dict)/len(ensg_ids):.1f}%)")
    print(f"  Saved to: {output_path}")
    
    # Also save as CSV for inspection
    mapping_df.to_csv(OUTPUT_DIR / "ensembl_to_entrez_geneformer.csv", index=False)
    
    # Check overlap with tAge features
    print("\nChecking overlap with tAge model...")
    import joblib
    model = joblib.load('/home/scroll/zzhang/transcriptome-aging/models/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl')
    # tAge features are numeric Entrez IDs stored in the pipeline
    # We know there are 10,487 features
    tage_features = set()  # We'll load actual features later
    
    # For now, just report the mapping coverage
    print(f"Entrez IDs mapped: {len(mapping_dict)}")


if __name__ == "__main__":
    main()
