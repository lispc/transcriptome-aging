"""
Evaluate best checkpoint on test set and save final model files.
"""

import json
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from datasets import load_from_disk
from transformers import BertForSequenceClassification, Trainer
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score
import shutil


def main():
    model_dir = Path("/home/scroll/zzhang/transcriptome-aging/models/geneformer_age_classifier_ddp")
    checkpoint_dir = model_dir / "checkpoint-4083"
    dataset_path = "/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_tokenized.dataset"

    # Load dataset
    print("Loading dataset...")
    ds = load_from_disk(dataset_path)
    label_list = sorted(set(ds["age_group"]))
    label2id = {l: i for i, l in enumerate(label_list)}
    id2label = {i: l for l, i in label2id.items()}

    def preprocess(example):
        example["labels"] = label2id[example["age_group"]]
        example["attention_mask"] = [1] * len(example["input_ids"])
        return example

    ds = ds.map(preprocess, num_proc=4)

    # Stratified split (same seed as training)
    donors = np.array(ds["donor_id"])
    labels = np.array(ds["labels"])
    donor_df = pd.DataFrame({"donor": donors, "label": labels})
    donor_labels = donor_df.groupby("donor")["label"].first().reset_index()

    _, temp_donors = train_test_split(
        donor_labels["donor"].values, test_size=0.3,
        stratify=donor_labels["label"].values, random_state=42,
    )
    val_donors, test_donors = train_test_split(
        temp_donors, test_size=0.5,
        stratify=donor_labels[donor_labels["donor"].isin(temp_donors)]["label"].values,
        random_state=42,
    )

    test_mask = np.isin(donors, test_donors)
    test_ds = ds.select(np.where(test_mask)[0])
    print(f"Test set: {len(test_ds)}")

    # Load model
    print("Loading model from best checkpoint...")
    model = BertForSequenceClassification.from_pretrained(str(checkpoint_dir))

    class Collator:
        def __init__(self, pad_token_id=0):
            self.pad_token_id = pad_token_id

        def __call__(self, features):
            max_len = max(len(f["input_ids"]) for f in features)
            batch = {"input_ids": [], "attention_mask": [], "labels": []}
            for f in features:
                ids = f["input_ids"]
                mask = f["attention_mask"]
                pad_len = max_len - len(ids)
                batch["input_ids"].append(ids + [self.pad_token_id] * pad_len)
                batch["attention_mask"].append(mask + [0] * pad_len)
                batch["labels"].append(f["labels"])
            return {k: torch.tensor(v) for k, v in batch.items()}

    def compute_metrics(eval_pred):
        logits, labels = eval_pred
        preds = np.argmax(logits, axis=1)
        return {
            "accuracy": accuracy_score(labels, preds),
            "f1_macro": f1_score(labels, preds, average="macro"),
            "f1_weighted": f1_score(labels, preds, average="weighted"),
        }

    trainer = Trainer(
        model=model,
        data_collator=Collator(),
        compute_metrics=compute_metrics,
    )

    print("Evaluating on test set...")
    results = trainer.evaluate(test_ds)
    print(f"Test results: {results}")

    # Save
    with open(model_dir / "test_metrics.json", "w") as f:
        json.dump(results, f, indent=2)
    with open(model_dir / "label_mapping.json", "w") as f:
        json.dump({"label2id": label2id, "id2label": id2label}, f, indent=2)

    # Copy best checkpoint files to root
    for fname in ["model.safetensors", "config.json"]:
        src = checkpoint_dir / fname
        dst = model_dir / fname
        shutil.copy2(src, dst)
        print(f"Copied {fname} to root")

    print("Done!")


if __name__ == "__main__":
    main()
