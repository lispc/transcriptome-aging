"""
In Silico Perturbation with fine-tuned Geneformer age classifier.
Supports checkpoint/resume: saves intermediate results every N genes.
"""

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


def predict_probs(model, cells, collator, device, batch_size=32):
    all_probs = []
    model.eval()
    with torch.no_grad():
        for i in range(0, len(cells), batch_size):
            batch = collator(cells[i:i + batch_size])
            batch = {k: v.to(device) for k, v in batch.items() if k != "labels"}
            outputs = model(**batch)
            probs = torch.softmax(outputs.logits, dim=-1).cpu().numpy()
            all_probs.append(probs)
    return np.vstack(all_probs)


def delete_gene(input_ids, gene_token_id):
    return [tid for tid in input_ids if tid != gene_token_id]


def main():
    model_dir = "/home/scroll/zzhang/transcriptome-aging/models/geneformer_age_classifier_ddp"
    dataset_path = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_tokenized.dataset"
    token_dict_path = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2/geneformer/token_dictionary_gc104M.pkl"
    top_genes_file = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/top2000_tokens.pkl"
    output_path = "/home/scroll/zzhang/transcriptome-aging/results/isp_finetuned_geneformer_age.csv"
    checkpoint_interval = 50  # save every N genes
    n_perturb_cells = 200
    batch_size = 32
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    print(f"Using device: {device}")

    # Load model
    print("Loading fine-tuned model...")
    model = BertForSequenceClassification.from_pretrained(model_dir)
    model.to(device)
    print(f"Model loaded: {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M params")

    # Load label mapping
    import json
    with open(Path(model_dir) / "label_mapping.json") as f:
        label_map = json.load(f)
    label2id = label_map["label2id"]
    id2label = {int(k): v for k, v in label_map["id2label"].items()}
    print(f"Labels: {id2label}")

    young_id = label2id.get("young", None)
    old_id = label2id.get("old", None)
    if young_id is None:
        raise ValueError("'young' not found in label mapping")

    # Load dataset
    print("Loading dataset...")
    ds = load_from_disk(dataset_path)

    # Sample cells for perturbation
    print(f"Sampling {n_perturb_cells} cells...")
    cells_per_group = n_perturb_cells // 3
    perturb_cells = []
    for group in ["young", "middle", "old"]:
        group_ds = ds.filter(lambda x: x["age_group"] == group, num_proc=4)
        n = min(cells_per_group, len(group_ds))
        sampled = group_ds.shuffle(seed=42).select(range(n))
        perturb_cells.extend([sampled[i] for i in range(n)])
    print(f"Perturbation cells: {len(perturb_cells)}")

    # Baseline predictions
    print("Computing baseline probabilities...")
    collator = GeneformerCollator(pad_token_id=0)
    baseline_probs = predict_probs(model, perturb_cells, collator, device, batch_size)
    print(f"Baseline young prob: {baseline_probs[:, young_id].mean():.3f}")
    if old_id is not None:
        print(f"Baseline old prob: {baseline_probs[:, old_id].mean():.3f}")

    # Load top genes
    with open(top_genes_file, "rb") as f:
        top_tokens = pickle.load(f)
    with open(token_dict_path, "rb") as f:
        token_dict = pickle.load(f)
    id2gene = {v: k for k, v in token_dict.items()}
    gene_tokens = [(tid, id2gene.get(tid, f"UNKNOWN_{tid}")) for tid in top_tokens if id2gene.get(tid, "").startswith("ENSG")]
    print(f"Genes to perturb: {len(gene_tokens)}")

    # Check for existing checkpoint
    checkpoint_path = Path(output_path).with_suffix(".checkpoint.csv")
    completed_genes = set()
    results = []

    if checkpoint_path.exists():
        print(f"Found checkpoint: {checkpoint_path}")
        df_ckpt = pd.read_csv(checkpoint_path)
        results = df_ckpt.to_dict("records")
        completed_genes = set(df_ckpt["ensembl_id"].tolist())
        print(f"Resuming from {len(completed_genes)} completed genes")

    # Filter to uncompleted genes
    remaining_genes = [(tid, ens) for tid, ens in gene_tokens if ens not in completed_genes]
    print(f"Remaining genes: {len(remaining_genes)}")

    # Run perturbations
    for idx, (tid, ens) in enumerate(tqdm(remaining_genes, desc="Perturbing genes")):
        n_present = sum(1 for c in perturb_cells if tid in c["input_ids"])
        if n_present == 0:
            results.append({
                "ensembl_id": ens, "token_id": tid,
                "delta_young_mean": 0, "delta_young_std": 0,
                "delta_old_mean": 0, "delta_old_std": 0,
                "n_cells_with_gene": 0,
            })
        else:
            perturbed = []
            for cell in perturb_cells:
                pcell = dict(cell)
                pcell["input_ids"] = delete_gene(cell["input_ids"], tid)
                pcell["attention_mask"] = [1] * len(pcell["input_ids"])
                perturbed.append(pcell)

            perturbed_probs = predict_probs(model, perturbed, collator, device, batch_size)
            delta_young = perturbed_probs[:, young_id] - baseline_probs[:, young_id]
            delta_old = perturbed_probs[:, old_id] - baseline_probs[:, old_id] if old_id is not None else np.zeros(len(delta_young))

            results.append({
                "ensembl_id": ens,
                "token_id": tid,
                "delta_young_mean": delta_young.mean(),
                "delta_young_std": delta_young.std(),
                "delta_old_mean": delta_old.mean(),
                "delta_old_std": delta_old.std(),
                "n_cells_with_gene": n_present,
            })

        # Save checkpoint periodically
        if (idx + 1) % checkpoint_interval == 0:
            pd.DataFrame(results).to_csv(checkpoint_path, index=False)
            print(f"  Checkpoint saved: {len(results)}/{len(gene_tokens)} genes")

    # Final save
    df = pd.DataFrame(results)
    df = df.sort_values("delta_young_mean", ascending=False)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    if checkpoint_path.exists():
        checkpoint_path.unlink()

    print(f"\nResults saved to {output_path}")
    print(f"\nTop 10 PRO-AGING genes (deletion -> younger; inhibit these):")
    print(df.head(10)[["ensembl_id", "delta_young_mean", "n_cells_with_gene"]].to_string())
    print(f"\nTop 10 ANTI-AGING genes (deletion -> older; activate these):")
    print(df.tail(10)[["ensembl_id", "delta_young_mean", "n_cells_with_gene"]].to_string())


if __name__ == "__main__":
    main()
