"""
In Silico Perturbation (ISP) for Geneformer age classifier.

For each gene in the Geneformer vocabulary, we perform "delete" perturbation
on a representative subset of cells and measure the change in predicted age class probabilities.

A gene that, when deleted, shifts predictions toward "young" is a potential pro-aging gene.
A gene that, when deleted, shifts predictions toward "old" is a potential anti-aging gene.
"""

import argparse
import json
import pickle
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from datasets import load_from_disk
from transformers import BertForSequenceClassification
from tqdm import tqdm


class GeneformerCollator:
    def __init__(self, pad_token_id=0, max_length=4096):
        self.pad_token_id = pad_token_id
        self.max_length = max_length

    def __call__(self, features):
        max_len = min(max(len(f["input_ids"]) for f in features), self.max_length)
        batch = {"input_ids": [], "attention_mask": [], "labels": []}
        for f in features:
            ids = f["input_ids"][:self.max_length]
            mask = [1] * len(ids)
            pad_len = max_len - len(ids)
            batch["input_ids"].append(ids + [self.pad_token_id] * pad_len)
            batch["attention_mask"].append(mask + [0] * pad_len)
            batch["labels"].append(f.get("labels", -1))
        return {k: torch.tensor(v) for k, v in batch.items()}


def delete_gene(input_ids, gene_token_id, pad_token_id=0):
    """Remove gene_token_id from input_ids sequence."""
    result = [tid for tid in input_ids if tid != gene_token_id]
    return result


def perturb_cells(cells, gene_token_id, collator, model, device, batch_size=32):
    """Run perturbation on a batch of cells and return probability changes."""
    perturbed = []
    for cell in cells:
        pcell = dict(cell)
        pcell["input_ids"] = delete_gene(cell["input_ids"], gene_token_id)
        pcell["attention_mask"] = [1] * len(pcell["input_ids"])
        perturbed.append(pcell)

    all_probs = []
    with torch.no_grad():
        for i in range(0, len(perturbed), batch_size):
            batch = collator(perturbed[i:i + batch_size])
            batch = {k: v.to(device) for k, v in batch.items()}
            outputs = model(**{k: v for k, v in batch.items() if k != "labels"})
            probs = torch.softmax(outputs.logits, dim=-1).cpu().numpy()
            all_probs.append(probs)

    return np.vstack(all_probs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", default="/home/scroll/zzhang/transcriptome-aging/models/geneformer_age_classifier")
    parser.add_argument("--dataset", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_tokenized.dataset")
    parser.add_argument("--token-dict", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2/geneformer/token_dictionary_gc104M.pkl")
    parser.add_argument("--output", default="/home/scroll/zzhang/transcriptome-aging/results/isp_geneformer_age.csv")
    parser.add_argument("--n-cells-per-group", type=int, default=100, help="Number of cells to perturb per age group")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load model
    print("Loading model...")
    model = BertForSequenceClassification.from_pretrained(args.model_dir)
    model.to(device)
    model.eval()

    # Load label mapping
    with open(Path(args.model_dir) / "label_mapping.json") as f:
        label_map = json.load(f)
    label2id = label_map["label2id"]
    id2label = label_map["id2label"]
    print(f"Labels: {id2label}")

    # Determine which label index corresponds to "young"
    young_id = label2id.get("young", None)
    old_id = label2id.get("old", None)
    if young_id is None:
        raise ValueError("'young' not found in label mapping")

    # Load dataset
    print("Loading dataset...")
    ds = load_from_disk(args.dataset)

    # Sample representative cells per age group
    cells_to_perturb = []
    for group in ["young", "middle", "old"]:
        group_ds = ds.filter(lambda x: x["age_group"] == group, num_proc=4)
        n = min(args.n_cells_per_group, len(group_ds))
        sampled = group_ds.shuffle(seed=42).select(range(n))
        cells_to_perturb.extend([sampled[i] for i in range(n)])
    print(f"Perturbing {len(cells_to_perturb)} cells")

    # Get baseline probabilities
    print("Computing baseline probabilities...")
    collator = GeneformerCollator(pad_token_id=0)
    baseline_probs = []
    with torch.no_grad():
        for i in tqdm(range(0, len(cells_to_perturb), args.batch_size)):
            batch = collator(cells_to_perturb[i:i + args.batch_size])
            batch = {k: v.to(device) for k, v in batch.items() if k != "labels"}
            outputs = model(**batch)
            probs = torch.softmax(outputs.logits, dim=-1).cpu().numpy()
            baseline_probs.append(probs)
    baseline_probs = np.vstack(baseline_probs)
    print(f"Baseline young prob mean: {baseline_probs[:, young_id].mean():.3f}")

    # Load token dictionary
    with open(args.token_dict, "rb") as f:
        token_dict = pickle.load(f)
    # Reverse: token_id -> ENSG
    id2gene = {v: k for k, v in token_dict.items()}

    # Filter to actual gene tokens (not special tokens)
    gene_tokens = [(tid, ens) for tid, ens in id2gene.items() if ens.startswith("ENSG")]
    print(f"Genes to perturb: {len(gene_tokens)}")

    # Run perturbations
    results = []
    for tid, ens in tqdm(gene_tokens, desc="Perturbing genes"):
        perturbed_probs = perturb_cells(
            cells_to_perturb, tid, collator, model, device, args.batch_size
        )

        # Impact metrics
        delta_young = perturbed_probs[:, young_id] - baseline_probs[:, young_id]
        delta_old = perturbed_probs[:, old_id] - baseline_probs[:, old_id] if old_id is not None else np.zeros(len(delta_young))

        results.append({
            "ensembl_id": ens,
            "token_id": tid,
            "delta_young_mean": delta_young.mean(),
            "delta_young_std": delta_young.std(),
            "delta_old_mean": delta_old.mean(),
            "delta_old_std": delta_old.std(),
            "n_cells_with_gene": sum(1 for c in cells_to_perturb if tid in c["input_ids"]),
        })

    # Save results
    df = pd.DataFrame(results)
    df = df.sort_values("delta_young_mean", ascending=False)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")
    print(f"\nTop 10 pro-aging genes (deletion -> younger):")
    print(df.head(10)[["ensembl_id", "delta_young_mean", "n_cells_with_gene"]].to_string())
    print(f"\nTop 10 anti-aging genes (deletion -> older):")
    print(df.tail(10)[["ensembl_id", "delta_young_mean", "n_cells_with_gene"]].to_string())


if __name__ == "__main__":
    main()
