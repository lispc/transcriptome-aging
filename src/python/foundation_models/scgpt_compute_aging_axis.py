"""
Compute aging axis with scGPT pretrained model (zero-shot).
Saves embeddings for later validation of Geneformer ISP results.
Uses GPU 1 to avoid conflict with ISP on GPU 0.

FIXES applied (2026-06-08):
- Set LD_LIBRARY_PATH so scGPT imports work (libstdc++ ABI issue).
- Correctly pass src_key_padding_mask to avoid padding-token contamination.
- Work around PyTorch 2.0+ nested-tensor bug on GPU by falling back to
  layer-by-layer forward when the full TransformerEncoder fails.
- Use scGPT's built-in _encode/_get_cell_emb_from_layer methods for
  consistency with pretraining.
"""

import os
# Fix scGPT import: miniforge has a newer libstdc++ than the system one.
os.environ["LD_LIBRARY_PATH"] = "/home/scroll/miniforge3/lib:" + os.environ.get("LD_LIBRARY_PATH", "")

import json
import numpy as np
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
    # Flash-attn stores QKV as fused Wqkv [3*embed_dim, embed_dim];
    # PyTorch MHA stores as fused in_proj [3*embed_dim, embed_dim].
    # Shapes match exactly, only the parameter name differs.
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
    print(f"Model loaded. Missing keys: {len(missing)}, Unexpected keys: {len(unexpected)}")
    if missing:
        print("  Missing (first 5):", missing[:5])
    if unexpected:
        print("  Unexpected (first 5):", unexpected[:5])

    # Verify attention weights actually loaded.
    ck_wqkv = checkpoint["transformer_encoder.layers.0.self_attn.in_proj_weight"]
    md_wqkv = dict(model.named_parameters())["transformer_encoder.layers.0.self_attn.in_proj_weight"]
    assert torch.allclose(ck_wqkv.cpu(), md_wqkv.cpu()), "Attention weight loading mismatch!"
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
    Encode a batch through the transformer, handling the PyTorch 2.0+
    nested-tensor bug on GPU.

    mask: bool tensor, True = real token, False = padding.
    """
    src = src.to(device)
    values = values.to(device)
    # PyTorch src_key_padding_mask convention: True = padding
    key_padding_mask = (~mask).to(device)

    gene_embs = model.encoder(src) + model.value_encoder(values)

    try:
        # Standard path (works on CPU; on GPU it may trigger nested-tensor optimisation)
        output = model.transformer_encoder(
            gene_embs, src_key_padding_mask=key_padding_mask
        )
    except (RuntimeError, TypeError) as e:
        err_msg = str(e).lower()
        if "nested" in err_msg or "to_padded_tensor" in err_msg:
            # Fallback: layer-by-layer avoids the nested-tensor code path.
            output = gene_embs
            for layer in model.transformer_encoder.layers:
                output = layer(output, src_key_padding_mask=key_padding_mask)
        else:
            raise

    # Cell embedding: use scGPT's CLS token (position 0).
    # mask is already bool, shape (batch, seq_len).
    # The CLS token is always real (position 0), but we still respect the mask.
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


def main():
    adata_path = "/home/scroll/zzhang/transcriptome-aging/data/scgpt/aida_v1_for_scgpt.h5ad"
    checkpoint_path = "/home/scroll/zzhang/transcriptome-aging/data/scgpt/weights/checkpoint-0623/scGPT_human/best_model.pt"
    vocab_path = "/home/scroll/zzhang/transcriptome-aging/data/scgpt/weights/checkpoint-0623/scGPT_human/vocab.json"
    output_dir = "/home/scroll/zzhang/transcriptome-aging/results"
    n_cells_per_group = 500
    batch_size = 32
    device = torch.device("cuda:1" if torch.cuda.device_count() > 1 else "cuda")

    print(f"Using device: {device}")

    print("Loading scGPT model...")
    model, vocab, model_args = load_scgpt_model(checkpoint_path, vocab_path, device)
    print(f"Model loaded: {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M params")

    print("Loading AIDA data...")
    adata = sc.read_h5ad(adata_path)
    print(f"Data shape: {adata.shape}")

    print(f"Sampling {n_cells_per_group} cells per age group...")
    young_adata = adata[adata.obs["age_group"] == "young", :]
    old_adata = adata[adata.obs["age_group"] == "old", :]

    n_young = min(n_cells_per_group, young_adata.n_obs)
    n_old = min(n_cells_per_group, old_adata.n_obs)

    np.random.seed(42)
    young_sample = young_adata[np.random.choice(young_adata.n_obs, n_young, replace=False), :]
    old_sample = old_adata[np.random.choice(old_adata.n_obs, n_old, replace=False), :]

    print(f"Young: {n_young}, Old: {n_old}")

    print("Extracting young embeddings...")
    young_embs = extract_embeddings(model, young_sample, vocab, device, batch_size=batch_size)
    print("Extracting old embeddings...")
    old_embs = extract_embeddings(model, old_sample, vocab, device, batch_size=batch_size)

    young_mean = young_embs.mean(axis=0)
    old_mean = old_embs.mean(axis=0)
    aging_axis = old_mean - young_mean
    aging_axis = aging_axis / np.linalg.norm(aging_axis)

    print(f"\n=== Results ===")
    print(f"Aging axis norm: {np.linalg.norm(aging_axis):.4f}")
    cos_sim = np.dot(old_mean / np.linalg.norm(old_mean), young_mean / np.linalg.norm(young_mean))
    print(f"Old-young cosine similarity: {cos_sim:.4f}")
    print(f"Old-young L2 distance: {np.linalg.norm(old_mean - young_mean):.4f}")

    # Per-dimension stats
    diffs = old_embs.mean(axis=0) - young_embs.mean(axis=0)
    print(f"Max |diff| across dims: {np.max(np.abs(diffs)):.4f}")
    print(f"Mean |diff| across dims: {np.mean(np.abs(diffs)):.4f}")
    print(f"Std of diff: {np.std(diffs):.4f}")

    # Save
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.save(out_dir / "scgpt_aging_axis.npy", aging_axis)
    np.save(out_dir / "scgpt_young_mean.npy", young_mean)
    np.save(out_dir / "scgpt_old_mean.npy", old_mean)
    np.save(out_dir / "scgpt_young_embs.npy", young_embs)
    np.save(out_dir / "scgpt_old_embs.npy", old_embs)
    print(f"\nSaved to {out_dir}")
    print("Done!")


if __name__ == "__main__":
    main()
