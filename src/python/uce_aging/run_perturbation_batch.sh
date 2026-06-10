#!/bin/bash
# Batch perturbation script - run on GPU 1 to avoid conflict with 33-layer on GPU 0
set -e

BASEDIR="/home/scroll/zzhang/transcriptome-aging"
MODEL="$BASEDIR/data/uce-repo/model_files/4layer_model.torch"
GENE_LIST="$BASEDIR/results/uce_spleen/perturbation_genes.txt"
WORKDIR="$BASEDIR/results/perturbation_work"
source ~/miniforge3/bin/activate aging-fm

cd "$BASEDIR"

# Prepare all perturbed files first
python -c "
import scanpy as sc
adata = sc.read_h5ad('results/uce_spleen/tms_spleen_subsampled_1k.h5ad')
if hasattr(adata.X, 'toarray'):
    adata.X = adata.X.toarray()
with open('$GENE_LIST') as f:
    genes = [l.strip() for l in f if l.strip()]
for gene in genes:
    adata_pert = adata.copy()
    if gene in adata_pert.var_names:
        idx = list(adata_pert.var_names).index(gene)
        adata_pert.X[:, idx] = 0
        adata_pert.write(f'$WORKDIR/perturbed_{gene}.h5ad')
        print(f'Prepared {gene}')
"

# Run UCE on each perturbed file
cd "$BASEDIR/data/uce-repo"
while IFS= read -r gene; do
    [ -z "$gene" ] && continue
    OUTDIR="$WORKDIR/uce_${gene}"
    mkdir -p "$OUTDIR"
    echo "=== Running UCE for $gene ==="
    CUDA_VISIBLE_DEVICES=1 python eval_single_anndata.py \
        --adata_path "$BASEDIR/results/perturbation_work/perturbed_${gene}.h5ad" \
        --dir "$OUTDIR" \
        --species mouse \
        --model_loc "$MODEL" \
        --batch_size 100 \
        --nlayers 4 \
        --output_dim 1280 2>&1 | tail -3
done < "$GENE_LIST"

echo "All perturbations complete!"
