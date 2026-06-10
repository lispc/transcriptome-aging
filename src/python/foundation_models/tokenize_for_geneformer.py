"""
Tokenize AIDA data for Geneformer V2-104M using TranscriptomeTokenizer.
"""

import sys
sys.path.insert(0, '/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2')

from geneformer import TranscriptomeTokenizer
from pathlib import Path


def main():
    input_file = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_for_geneformer.h5ad"
    output_dir = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized"
    output_prefix = "aida_v1_tokenized"
    
    # Custom attributes to preserve
    custom_attrs = {
        "age_group": "age_group",
        "age_numeric": "age_numeric", 
        "cell_type": "cell_type",
        "donor_id": "donor_id",
    }
    
    print("Initializing TranscriptomeTokenizer...")
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict=custom_attrs,
        nproc=4,
        model_version="V2",
        use_h5ad_index=True,
        
    )
    
    print(f"Tokenizing {input_file}...")
    tk.tokenize_data(
        data_directory=Path(input_file).parent,
        output_directory=output_dir,
        output_prefix=output_prefix,
        file_format="h5ad",
        input_identifier="aida_v1_for_geneformer",
    )
    
    print(f"Tokenized dataset saved to {output_dir}/{output_prefix}.dataset")


if __name__ == "__main__":
    main()
