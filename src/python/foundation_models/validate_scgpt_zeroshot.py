"""
Zero-shot validation with scGPT pretrained model.

For top candidate genes from Geneformer ISP, we check if their deletion
also shifts cells toward "young" in scGPT's embedding space.

FIXES applied (2026-06-08):
- Set LD_LIBRARY_PATH so scGPT imports work.
- Convert FlashMHA Wqkv -> standard MultiheadAttention in_proj weights.
- Use use_fast_transformer=False since flash-attn is not installed.
- Avoid nested-tensor bug on GPU by using layer-by-layer fallback.
- Extract correct cell embeddings from transformer output (not decoder predictions).
"""

import os
os.environ["LD_LIBRARY_PATH"] = "/home/scroll/miniforge3/lib:" + os.environ.get("LD_LIBRARY_PATH", "")

import argparse
import json
import numpy as np
import pandas as pd
import torch
import scanpy as sc
from pathlib import Path
from tqdm import tqdm


def load_scgpt_model(checkpoint_path, vocab_path, device):
    import sys
    sys.path.insert(0, '/home/scroll/zzhang/transcriptome-aging/data/scgpt/weights/scgpt-repo')
    from scgpt.model import TransformerModel
    from scgpt.tokenizer import GeneVocab

    vocab = GeneVocab.from_file(vocab_path)
    vocab.set_default_index(vocab["<pad>"])

    checkpoint = torch.load(checkpoint_path, map_location=device)

    # Convert FlashMHA Wqkv weights to standard MultiheadAttention in_proj weights.
    converted = {}
    for k, v in checkpoint.items():
        if "self_attn.Wqkv.weight" in k:
            new_k = k.replace("self_attn.Wqkv.weight", "self_attn.in_proj_weight")
            converted[new_k] = v
        elif "self_attn.Wqkv.bias" in k:
            new_k = k.replace("self_attn.Wqkv.bias", "self_attn.in_proj_bias")
            converted[new_k] = v
        else:
            converted[k] = v
    checkpoint = converted

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
        use_fast_transformer=False,  # flash-attn not available
        fast_transformer_backend="flash",
        pre_norm=False,
    )
    missing, unexpected = model.load_state_dict(checkpoint, strict=False)
    print(f"Model loaded. Missing: {len(missing)}, Unexpected: {len(unexpected)}")

    # Sanity check attention weights
    ck = checkpoint["transformer_encoder.layers.0.self_attn.in_proj_weight"]
    md = dict(model.named_parameters())["transformer_encoder.layers.0.self_attn.in_proj_weight"]
    assert torch.allclose(ck.cpu(), md.cpu()), "Attention weight loading mismatch!"
    print("Attention weight sanity check: PASSED")

    model.to(device)
    model.eval()
    return model, vocab, args


def prepare_cell_for_scgpt(cell_expr, cell_genes, vocab, max_seq_len=1200):
    sorted_idx = np.argsort(-cell_expr)
    sorted_genes = [cell_genes[i] for i in sorted_idx]
    sorted_values = cell_expr[sorted_idx]

    gene_ids = []
    values = []
    for g, v in zip(sorted_genes, sorted_values):
        if g in vocab:
            gene_ids.append(vocab[g])
            values.append(v)

    gene_ids = gene_ids[:max_seq_len]
    values = values[:max_seq_len]
    pad_len = max_seq_len - len(gene_ids)
    gene_ids = gene_ids + [vocab["<pad>"]] * pad_len
    values = values + [-2] * pad_len
    mask = [1] * (max_seq_len - pad_len) + [0] * pad_len

    return (
        torch.tensor(gene_ids, dtype=torch.long),
        torch.tensor(values, dtype=torch.float32),
        torch.tensor(mask, dtype=torch.bool),
    )


def safe_encode(model, src, values, mask, device):
    """
    Encode a batch through the transformer, handling nested-tensor bug on GPU.
    mask: bool tensor, True = real token, False = padding.
    Returns cell embeddings using CLS style (position 0).
    """
    src = src.to(device)
    values = values.to(device)
    key_padding_mask = (~mask).to(device)  # True = padding

    gene_embs = model.encoder(src) + model.value_encoder(values)

    try:
        output = model.transformer_encoder(
            gene_embs, src_key_padding_mask=key_padding_mask
        )
    except (RuntimeError, TypeError) as e:
        err_msg = str(e).lower()
        if "nested" in err_msg or "to_padded_tensor" in err_msg:
            output = gene_embs
            for layer in model.transformer_encoder.layers:
                output = layer(output, src_key_padding_mask=key_padding_mask)
        else:
            raise

    cell_emb = model._get_cell_emb_from_layer(output, weights=None)
    return cell_emb


def extract_embeddings(model, adata, vocab, device, max_seq_len=1200, batch_size=32):
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

            src = torch.stack(batch_src)
            vals = torch.stack(batch_vals)
            mask = torch.stack(batch_mask)

            cell_emb = safe_encode(model, src, vals, mask, device)
            embs.append(cell_emb.cpu().numpy())

    return np.vstack(embs)


def delete_gene_and_embed(cell_expr, cell_genes, gene_to_delete, model, vocab, device, max_seq_len=1200):
    """Delete a gene from a cell and extract embedding."""
    expr = cell_expr.copy()
    if gene_to_delete in cell_genes:
        idx = cell_genes.index(gene_to_delete)
        expr[idx] = 0

    src, vals, mask = prepare_cell_for_scgpt(expr, cell_genes, vocab, max_seq_len)
    src = src.unsqueeze(0)
    vals = vals.unsqueeze(0)
    mask = mask.unsqueeze(0)

    with torch.no_grad():
        cell_emb = safe_encode(model, src, vals, mask, device)

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

    print("Loading scGPT model...")
    model, vocab, model_args = load_scgpt_model(args.checkpoint, args.vocab, device)
    print(f"Model loaded: {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M params")
    print(f"Vocab size: {len(vocab)}")

    print("Loading AIDA data...")
    adata = sc.read_h5ad(args.adata)
    print(f"Data shape: {adata.shape}")

    print(f"Sampling {args.n_cells_per_group} cells per age group...")
    young_adata = adata[adata.obs["age_group"] == "young", :]
    old_adata = adata[adata.obs["age_group"] == "old", :]

    n_young = min(args.n_cells_per_group, young_adata.n_obs)
    n_old = min(args.n_cells_per_group, old_adata.n_obs)

    np.random.seed(42)
    young_sample = young_adata[np.random.choice(young_adata.n_obs, n_young, replace=False), :]
    old_sample = old_adata[np.random.choice(old_adata.n_obs, n_old, replace=False), :]

    print("Extracting young embeddings...")
    young_embs = extract_embeddings(model, young_sample, vocab, device, args.batch_size)
    print("Extracting old embeddings...")
    old_embs = extract_embeddings(model, old_sample, vocab, device, args.batch_size)

    young_mean = young_embs.mean(axis=0)
    old_mean = old_embs.mean(axis=0)
    aging_axis = old_mean - young_mean
    aging_axis = aging_axis / np.linalg.norm(aging_axis)
    print(f"Aging axis norm: {np.linalg.norm(aging_axis):.4f}")
    cos_sim = np.dot(old_mean / np.linalg.norm(old_mean), young_mean / np.linalg.norm(young_mean))
    print(f"Old-young cosine similarity: {cos_sim:.4f}")

    print("Loading Geneformer ISP results...")
    isp_df = pd.read_csv(args.isp_results)
    top_genes = isp_df.head(args.top_n_genes)["ensembl_id"].tolist()
    print(f"Top {len(top_genes)} genes to validate: {top_genes[:5]}...")

    with open("/home/scroll/zzhang/transcriptome-aging/data/gene_id_mappings/ensg_to_symbol_scgpt.json") as f:
        ensg_to_symbol = json.load(f)

    n_perturb = 100
    perturb_adata = adata[np.random.choice(adata.n_obs, min(n_perturb, adata.n_obs), replace=False), :]
    gene_names = list(perturb_adata.var.index)

    print("Computing baseline embeddings...")
    baseline_embs = extract_embeddings(model, perturb_adata, vocab, device, args.batch_size)
    baseline_proj = np.dot(baseline_embs, aging_axis)

    results = []
    for ens in tqdm(top_genes, desc="Validating genes"):
        symbol = ensg_to_symbol.get(ens, None)
        if symbol is None or symbol not in vocab:
            results.append({
                "ensembl_id": ens, "gene_symbol": symbol or "N/A",
                "in_scgpt_vocab": False, "delta_proj_mean": np.nan,
                "delta_proj_std": np.nan,
                "n_cells_with_gene": np.nan,
            })
            continue

        perturbed_embs = []
        n_cells_with_gene = 0
        for i in range(perturb_adata.n_obs):
            expr = perturb_adata.X[i, :]
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()
            else:
                expr = np.array(expr).flatten()

            if symbol in gene_names and expr[gene_names.index(symbol)] > 0:
                n_cells_with_gene += 1

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
            "n_cells_with_gene": n_cells_with_gene,
        })

    df = pd.DataFrame(results)
    df = df.sort_values("delta_proj_mean", ascending=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    print(f"\nResults saved to {args.output}")
    print(f"\nTop validated pro-aging genes (scGPT agrees with Geneformer, deletion -> younger):")
    print(df.head(10)[["ensembl_id", "gene_symbol", "delta_proj_mean", "delta_proj_std"]].to_string())
    print(f"\nGenes where scGPT disagrees with Geneformer (deletion -> older):")
    print(df.tail(5)[["ensembl_id", "gene_symbol", "delta_proj_mean", "delta_proj_std"]].to_string())


if __name__ == "__main__":
    main()
