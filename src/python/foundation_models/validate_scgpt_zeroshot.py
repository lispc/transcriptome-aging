"""
Zero-shot validation with scGPT pretrained model.

For top candidate genes from Geneformer ISP, we check if their deletion
also shifts cells toward "young" in scGPT's embedding space.

Strategy:
1. Load pretrained scGPT model
2. Convert AIDA data to scGPT format (top 1200 genes + continuous values)
3. Compute aging axis = mean(old_emb) - mean(young_emb)
4. For each candidate gene, delete it and re-extract cell embedding
5. Measure projection shift onto aging axis
"""

import argparse
import json
import pickle
import numpy as np
import pandas as pd
import torch
import scanpy as sc
from pathlib import Path
from tqdm import tqdm


def load_scgpt_model(checkpoint_path, vocab_path, device):
    """Load pretrained scGPT model from checkpoint."""
    import sys
    sys.path.insert(0, '/home/scroll/zzhang/transcriptome-aging/data/scgpt/weights/scgpt-repo')
    from scgpt.model import TransformerModel
    from scgpt.tokenizer import GeneVocab

    # Load vocab
    vocab = GeneVocab.from_file(vocab_path)
    vocab.set_default_index(vocab["<pad>"])

    # Load checkpoint
    checkpoint = torch.load(checkpoint_path, map_location=device)

    # Model args from checkpoint or use defaults from args.json
    with open(Path(checkpoint_path).parent / "args.json") as f:
        args = json.load(f)

    model = TransformerModel(
        ntoken=len(vocab),
        d_model=args.get("embsize", 512),
        nhead=args.get("nheads", 8),
        d_hid=args.get("d_hid", 512),
        nlayers=args.get("nlayers", 12),
        nlayers_cls=args.get("n_layers_cls", 3),
        n_cls=1,
        vocab=vocab,
        dropout=args.get("dropout", 0.2),
        pad_token="<pad>",
        pad_value=-2,
        do_mvc=False,
        do_dab=False,
        use_batch_labels=False,
        input_emb_style="continuous",
        n_input_bins=51,
        cell_emb_style="cls",
        mvc_decoder_style="inner product",
        ecs_threshold=0.3,
        explicit_zero_prob=False,
        use_fast_transformer=args.get("fast_transformer", True),
        fast_transformer_backend="flash",
        pre_norm=False,
    )

    # Load weights
    model.load_state_dict(checkpoint)
    model.to(device)
    model.eval()

    return model, vocab, args


def prepare_cell_for_scgpt(cell_expr, cell_genes, vocab, max_seq_len=1200):
    """
    Convert a single cell's expression to scGPT input format.
    cell_expr: array of expression values
    cell_genes: list of gene symbols
    """
    # Sort genes by expression (descending)
    sorted_idx = np.argsort(-cell_expr)
    sorted_genes = [cell_genes[i] for i in sorted_idx]
    sorted_values = cell_expr[sorted_idx]

    # Filter to genes in vocab
    gene_ids = []
    values = []
    for g, v in zip(sorted_genes, sorted_values):
        if g in vocab:
            gene_ids.append(vocab[g])
            values.append(v)

    # Truncate to max_seq_len
    gene_ids = gene_ids[:max_seq_len]
    values = values[:max_seq_len]

    # Pad
    pad_len = max_seq_len - len(gene_ids)
    gene_ids = gene_ids + [vocab["<pad>"]] * pad_len
    values = values + [-2] * pad_len  # pad_value = -2
    mask = [1] * (max_seq_len - pad_len) + [0] * pad_len

    return torch.tensor(gene_ids), torch.tensor(values, dtype=torch.float32), torch.tensor(mask)


def extract_embeddings(model, adata, vocab, device, max_seq_len=1200, batch_size=32):
    """Extract cell embeddings from scGPT."""
    gene_names = list(adata.var.index)
    embs = []

    model.eval()
    with torch.no_grad():
        for start in tqdm(range(0, adata.n_obs, batch_size), desc="Extracting embeddings"):
            end = min(start + batch_size, adata.n_obs)

            batch_src = []
            batch_vals = []
            batch_mask = []

            for i in range(start, end):
                expr = adata.X[i, :]
                if hasattr(expr, "toarray"):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.array(expr).flatten()

                src, vals, mask = prepare_cell_for_scgpt(expr, gene_names, vocab, max_seq_len)
                batch_src.append(src)
                batch_vals.append(vals)
                batch_mask.append(mask)

            src = torch.stack(batch_src).to(device)
            vals = torch.stack(batch_vals).to(device)
            mask = torch.stack(batch_mask).to(device)

            # scGPT forward returns a dict with 'cell_emb' when CLS=True
            # But since no_cls=true in training, we use average pooling instead
            outputs = model(src, vals, src_key_padding_mask=~mask.bool())

            # Get cell embedding: average of gene embeddings (excluding padding)
            gene_embs = outputs["mlm_output"] if "mlm_output" in outputs else outputs[0]
            # gene_embs: [batch, seq_len, hidden]

            # Mean pool over non-padding positions
            mask_expanded = mask.unsqueeze(-1).float()
            cell_emb = (gene_embs * mask_expanded).sum(dim=1) / mask_expanded.sum(dim=1)
            embs.append(cell_emb.cpu().numpy())

    return np.vstack(embs)


def delete_gene_and_embed(cell_expr, cell_genes, gene_to_delete, model, vocab, device, max_seq_len=1200):
    """Delete a gene from a cell and extract embedding."""
    # Zero out the gene's expression
    expr = cell_expr.copy()
    if gene_to_delete in cell_genes:
        idx = cell_genes.index(gene_to_delete)
        expr[idx] = 0

    src, vals, mask = prepare_cell_for_scgpt(expr, cell_genes, vocab, max_seq_len)
    src = src.unsqueeze(0).to(device)
    vals = vals.unsqueeze(0).to(device)
    mask = mask.unsqueeze(0).to(device)

    with torch.no_grad():
        outputs = model(src, vals, src_key_padding_mask=~mask.bool())
        gene_embs = outputs["mlm_output"] if "mlm_output" in outputs else outputs[0]
        mask_expanded = mask.unsqueeze(-1).float()
        cell_emb = (gene_embs * mask_expanded).sum(dim=1) / mask_expanded.sum(dim=1)

    return cell_emb.cpu().numpy()[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", default="/home/scroll/zzhang/transcriptome-aging/data/scgpt/aida_v1_for_scgpt.h5ad")
    parser.add_argument("--checkpoint", default="/home/scroll/zzhang/transcriptome-aging/data/scgpt/weights/checkpoint-0623/scGPT_human/best_model.pt")
    parser.add_argument("--vocab", default="/home/scroll/zzhang/transcriptome-aging/data/scgpt/weights/checkpoint-0623/scGPT_human/vocab.json")
    parser.add_argument("--isp-results", default="/home/scroll/zzhang/transcriptome-aging/results/isp_finetuned_geneformer_age.csv")
    parser.add_argument("--output", default="/home/scroll/zzhang/transcriptome-aging/results/scgpt_validation.csv")
    parser.add_argument("--n-cells-per-group", type=int, default=200)
    parser.add_argument("--top-n-genes", type=int, default=50)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load scGPT model
    print("Loading scGPT model...")
    model, vocab, model_args = load_scgpt_model(args.checkpoint, args.vocab, device)
    print(f"Model loaded: {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M params")
    print(f"Vocab size: {len(vocab)}")

    # Load AIDA data
    print("Loading AIDA data...")
    adata = sc.read_h5ad(args.adata)
    print(f"Data shape: {adata.shape}")

    # Sample young/old cells for aging axis
    print(f"Sampling {args.n_cells_per_group} cells per age group...")
    young_adata = adata[adata.obs["age_group"] == "young", :]
    old_adata = adata[adata.obs["age_group"] == "old", :]

    n_young = min(args.n_cells_per_group, young_adata.n_obs)
    n_old = min(args.n_cells_per_group, old_adata.n_obs)

    young_sample = young_adata[np.random.choice(young_adata.n_obs, n_young, replace=False), :]
    old_sample = old_adata[np.random.choice(old_adata.n_obs, n_old, replace=False), :]

    # Extract embeddings
    print("Extracting young embeddings...")
    young_embs = extract_embeddings(model, young_sample, vocab, device, args.batch_size)
    print("Extracting old embeddings...")
    old_embs = extract_embeddings(model, old_sample, vocab, device, args.batch_size)

    young_mean = young_embs.mean(axis=0)
    old_mean = old_embs.mean(axis=0)
    aging_axis = old_mean - young_mean
    aging_axis = aging_axis / np.linalg.norm(aging_axis)
    print(f"Aging axis norm: {np.linalg.norm(aging_axis):.4f}")

    # Load Geneformer ISP results
    print("Loading Geneformer ISP results...")
    isp_df = pd.read_csv(args.isp_results)
    # Top pro-aging genes (deletion -> younger)
    top_genes = isp_df.head(args.top_n_genes)["ensembl_id"].tolist()
    print(f"Top {len(top_genes)} genes to validate: {top_genes[:5]}...")

    # Map ENSG to gene symbols
    with open("/home/scroll/zzhang/transcriptome-aging/data/gene_id_mappings/ensg_to_symbol_scgpt.json") as f:
        ensg_to_symbol = json.load(f)

    # Sample cells for perturbation
    n_perturb = 100
    perturb_adata = adata[np.random.choice(adata.n_obs, min(n_perturb, adata.n_obs), replace=False), :]
    gene_names = list(perturb_adata.var.index)

    # Baseline embeddings
    print("Computing baseline embeddings...")
    baseline_embs = extract_embeddings(model, perturb_adata, vocab, device, args.batch_size)
    baseline_proj = np.dot(baseline_embs, aging_axis)

    # Validate each candidate gene
    results = []
    for ens in tqdm(top_genes, desc="Validating genes"):
        symbol = ensg_to_symbol.get(ens, None)
        if symbol is None or symbol not in vocab:
            results.append({
                "ensembl_id": ens, "gene_symbol": symbol or "N/A",
                "in_scgpt_vocab": False, "delta_proj": np.nan,
            })
            continue

        perturbed_embs = []
        for i in range(perturb_adata.n_obs):
            expr = perturb_adata.X[i, :]
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()
            else:
                expr = np.array(expr).flatten()

            emb = delete_gene_and_embed(expr, gene_names, symbol, model, vocab, device)
            perturbed_embs.append(emb)

        perturbed_embs = np.vstack(perturbed_embs)
        perturbed_proj = np.dot(perturbed_embs, aging_axis)
        delta = perturbed_proj - baseline_proj

        results.append({
            "ensembl_id": ens,
            "gene_symbol": symbol,
            "in_scgpt_vocab": True,
            "delta_proj_mean": delta.mean(),
            "delta_proj_std": delta.std(),
            "n_cells_with_gene": sum(1 for i in range(perturb_adata.n_obs) if symbol in gene_names and perturb_adata.X[i, gene_names.index(symbol)] > 0),
        })

    # Save
    df = pd.DataFrame(results)
    df = df.sort_values("delta_proj_mean", ascending=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    print(f"\nResults saved to {args.output}")
    print(f"\nTop validated pro-aging genes (scGPT agrees with Geneformer):")
    print(df.head(10)[["ensembl_id", "gene_symbol", "delta_proj_mean"]].to_string())
    print(f"\nGenes where scGPT disagrees with Geneformer:")
    print(df.tail(5)[["ensembl_id", "gene_symbol", "delta_proj_mean"]].to_string())


if __name__ == "__main__":
    main()
