"""
Zero-shot In Silico Perturbation with Geneformer (no fine-tuning required).

Strategy:
1. Load pretrained Geneformer (BertModel) to extract CLS embeddings
2. Compute "aging axis" as mean(CLS_old) - mean(CLS_young)
3. For each gene, delete it from cell's token sequence, re-extract CLS
4. Measure projection of CLS shift onto aging axis
5. Genes whose deletion shifts cells toward "young" are pro-aging
   Genes whose deletion shifts cells toward "old" are anti-aging
"""

import argparse
import pickle
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from datasets import load_from_disk
from transformers import BertModel
from tqdm import tqdm


class GeneformerCollator:
    def __init__(self, pad_token_id=0, max_length=4096):
        self.pad_token_id = pad_token_id
        self.max_length = max_length

    def __call__(self, features):
        max_len = min(max(len(f["input_ids"]) for f in features), self.max_length)
        batch = {"input_ids": [], "attention_mask": []}
        for f in features:
            ids = f["input_ids"][:self.max_length]
            mask = [1] * len(ids)
            pad_len = max_len - len(ids)
            batch["input_ids"].append(ids + [self.pad_token_id] * pad_len)
            batch["attention_mask"].append(mask + [0] * pad_len)
        return {k: torch.tensor(v) for k, v in batch.items()}


def extract_cls(model, cells, collator, device, batch_size=128):
    """Extract CLS embeddings for a list of cells."""
    embs = []
    model.eval()
    with torch.no_grad():
        for i in range(0, len(cells), batch_size):
            batch = collator(cells[i:i + batch_size])
            batch = {k: v.to(device) for k, v in batch.items()}
            outputs = model(**batch)
            cls_emb = outputs.last_hidden_state[:, 0, :].cpu().numpy()
            embs.append(cls_emb)
    return np.vstack(embs)


def delete_gene(input_ids, gene_token_id):
    return [tid for tid in input_ids if tid != gene_token_id]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-path", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2/Geneformer-V2-104M")
    parser.add_argument("--dataset", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_tokenized.dataset")
    parser.add_argument("--token-dict", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2/geneformer/token_dictionary_gc104M.pkl")
    parser.add_argument("--top-genes-file", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/top2000_tokens.pkl")
    parser.add_argument("--output", default="/home/scroll/zzhang/transcriptome-aging/results/isp_zeroshot_geneformer.csv")
    parser.add_argument("--n-cells-per-group", type=int, default=200)
    parser.add_argument("--n-perturb-cells", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    print("Loading pretrained Geneformer...")
    model = BertModel.from_pretrained(args.model_path)
    model.to(device)
    model.eval()
    print(f"Model loaded: {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M params")

    print("Loading dataset...")
    ds = load_from_disk(args.dataset)

    # Sample cells for aging axis
    print(f"Sampling {args.n_cells_per_group} cells per age group...")
    young_cells = []
    old_cells = []
    for group, target_list in [("young", young_cells), ("old", old_cells)]:
        group_ds = ds.filter(lambda x: x["age_group"] == group, num_proc=4)
        n = min(args.n_cells_per_group, len(group_ds))
        sampled = group_ds.shuffle(seed=42).select(range(n))
        target_list.extend([sampled[i] for i in range(n)])
    print(f"Young: {len(young_cells)}, Old: {len(old_cells)}")

    collator = GeneformerCollator(pad_token_id=0)

    # Compute aging axis
    print("Extracting CLS embeddings for aging axis...")
    young_embs = extract_cls(model, young_cells, collator, device, args.batch_size)
    old_embs = extract_cls(model, old_cells, collator, device, args.batch_size)

    young_mean = young_embs.mean(axis=0)
    old_mean = old_embs.mean(axis=0)
    aging_axis = old_mean - young_mean
    aging_axis = aging_axis / np.linalg.norm(aging_axis)
    print(f"Aging axis norm: {np.linalg.norm(aging_axis):.4f}")
    cos_sim = np.dot(old_mean / np.linalg.norm(old_mean), young_mean / np.linalg.norm(young_mean))
    print(f"Old-young cosine similarity: {cos_sim:.4f}")

    # Sample cells for perturbation
    print("Sampling cells for perturbation...")
    perturb_ds = ds.shuffle(seed=42).select(range(min(args.n_perturb_cells, len(ds))))
    perturb_cells = [perturb_ds[i] for i in range(len(perturb_ds))]
    print(f"Perturbation cells: {len(perturb_cells)}")

    # Baseline embeddings
    print("Computing baseline embeddings...")
    baseline_embs = extract_cls(model, perturb_cells, collator, device, args.batch_size)
    baseline_proj = np.dot(baseline_embs, aging_axis)
    print(f"Baseline projection range: {baseline_proj.min():.3f} to {baseline_proj.max():.3f}")

    # Load top genes to perturb
    with open(args.top_genes_file, "rb") as f:
        top_tokens = pickle.load(f)
    with open(args.token_dict, "rb") as f:
        token_dict = pickle.load(f)
    id2gene = {v: k for k, v in token_dict.items()}
    gene_tokens = [(tid, id2gene.get(tid, f"UNKNOWN_{tid}")) for tid in top_tokens if id2gene.get(tid, "").startswith("ENSG")]
    print(f"Genes to perturb: {len(gene_tokens)}")

    # Run perturbations
    results = []
    for tid, ens in tqdm(gene_tokens, desc="Perturbing genes"):
        n_present = sum(1 for c in perturb_cells if tid in c["input_ids"])
        if n_present == 0:
            results.append({
                "ensembl_id": ens, "token_id": tid,
                "delta_proj_mean": 0, "delta_proj_std": 0,
                "n_cells_with_gene": 0,
            })
            continue

        perturbed = []
        for cell in perturb_cells:
            pcell = dict(cell)
            pcell["input_ids"] = delete_gene(cell["input_ids"], tid)
            pcell["attention_mask"] = [1] * len(pcell["input_ids"])
            perturbed.append(pcell)

        perturbed_embs = extract_cls(model, perturbed, collator, device, args.batch_size)
        perturbed_proj = np.dot(perturbed_embs, aging_axis)
        delta = perturbed_proj - baseline_proj

        results.append({
            "ensembl_id": ens,
            "token_id": tid,
            "delta_proj_mean": delta.mean(),
            "delta_proj_std": delta.std(),
            "n_cells_with_gene": n_present,
        })

    # Save and rank
    df = pd.DataFrame(results)
    df = df.sort_values("delta_proj_mean", ascending=True)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    print(f"\nResults saved to {args.output}")
    print(f"\nTop 10 PRO-AGING genes (deletion -> younger):")
    print(df.head(10)[["ensembl_id", "delta_proj_mean", "n_cells_with_gene"]].to_string())
    print(f"\nTop 10 ANTI-AGING genes (deletion -> older):")
    print(df.tail(10)[["ensembl_id", "delta_proj_mean", "n_cells_with_gene"]].to_string())

    # Save aging axis
    out_dir = Path(args.output).parent
    np.save(out_dir / "aging_axis.npy", aging_axis)
    np.save(out_dir / "young_mean.npy", young_mean)
    np.save(out_dir / "old_mean.npy", old_mean)
    print("\nAging axis vectors saved.")


if __name__ == "__main__":
    main()
