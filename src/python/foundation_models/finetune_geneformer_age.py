"""
Fine-tune Geneformer V2-104M for age group classification (young/middle/old).

Strategy:
- Load pre-trained Geneformer as BertForSequenceClassification
- Freeze first N layers for efficient fine-tuning
- Stratified split by donor_id to avoid data leakage
- Train with HuggingFace Trainer API
"""

import argparse
import json
import pickle
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from datasets import load_from_disk
from sklearn.model_selection import StratifiedGroupKFold
from transformers import (
    BertForSequenceClassification,
    TrainingArguments,
    Trainer,
    EarlyStoppingCallback,
)
from transformers import DataCollatorWithPadding


def freeze_layers(model, n_freeze):
    """Freeze first n_freeze encoder layers."""
    for param in model.bert.embeddings.parameters():
        param.requires_grad = False
    for i, layer in enumerate(model.bert.encoder.layer):
        if i < n_freeze:
            for param in layer.parameters():
                param.requires_grad = False
    print(f"Frozen first {n_freeze} encoder layers + embeddings")


def compute_metrics(eval_pred):
    from sklearn.metrics import accuracy_score, f1_score
    logits, labels = eval_pred
    preds = np.argmax(logits, axis=1)
    return {
        "accuracy": accuracy_score(labels, preds),
        "f1_macro": f1_score(labels, preds, average="macro"),
        "f1_weighted": f1_score(labels, preds, average="weighted"),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_tokenized.dataset")
    parser.add_argument("--model-path", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2/Geneformer-V2-104M")
    parser.add_argument("--output-dir", default="/home/scroll/zzhang/transcriptome-aging/models/geneformer_age_classifier")
    parser.add_argument("--freeze-layers", type=int, default=6, help="Freeze first N transformer layers")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=5e-5)
    parser.add_argument("--weight-decay", type=float, default=0.01)
    parser.add_argument("--warmup-ratio", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # Load dataset
    print("Loading dataset...")
    ds = load_from_disk(args.dataset)
    print(f"Dataset size: {len(ds)}")

    # Encode labels
    label_list = sorted(set(ds["age_group"]))
    label2id = {l: i for i, l in enumerate(label_list)}
    id2label = {i: l for l, i in label2id.items()}
    print(f"Labels: {label_list}")

    # Add labels and attention_mask
    def preprocess(example):
        example["labels"] = label2id[example["age_group"]]
        # Create attention_mask: 1 for real tokens, 0 for padding
        # Since there's no padding yet, all are 1
        example["attention_mask"] = [1] * len(example["input_ids"])
        return example

    ds = ds.map(preprocess, num_proc=4)

    # Stratified split by donor_id
    donors = np.array(ds["donor_id"])
    labels = np.array(ds["labels"])
    
    # Use donor-level stratification: each donor gets one representative label
    donor_df = pd.DataFrame({"donor": donors, "label": labels})
    donor_labels = donor_df.groupby("donor")["label"].first().reset_index()
    
    from sklearn.model_selection import train_test_split
    train_donors, temp_donors = train_test_split(
        donor_labels["donor"].values,
        test_size=0.3,
        stratify=donor_labels["label"].values,
        random_state=args.seed,
    )
    val_donors, test_donors = train_test_split(
        temp_donors,
        test_size=0.5,
        stratify=donor_labels[donor_labels["donor"].isin(temp_donors)]["label"].values,
        random_state=args.seed,
    )

    train_mask = np.isin(donors, train_donors)
    val_mask = np.isin(donors, val_donors)
    test_mask = np.isin(donors, test_donors)

    train_ds = ds.select(np.where(train_mask)[0])
    val_ds = ds.select(np.where(val_mask)[0])
    test_ds = ds.select(np.where(test_mask)[0])

    print(f"Train: {len(train_ds)}, Val: {len(val_ds)}, Test: {len(test_ds)}")
    print(f"Train donors: {len(train_donors)}, Val donors: {len(val_donors)}, Test donors: {len(test_donors)}")

    # Load model
    print("Loading model...")
    model = BertForSequenceClassification.from_pretrained(
        args.model_path,
        num_labels=len(label_list),
        id2label=id2label,
        label2id=label2id,
    )

    if args.freeze_layers > 0:
        freeze_layers(model, args.freeze_layers)

    # Count trainable parameters
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    total = sum(p.numel() for p in model.parameters())
    print(f"Trainable params: {trainable/1e6:.1f}M / {total/1e6:.1f}M")

    # Training arguments
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    training_args = TrainingArguments(
        output_dir=str(output_dir),
        num_train_epochs=3,
        per_device_train_batch_size=8,
        per_device_eval_batch_size=16,
        gradient_accumulation_steps=4,
        fp16=True,
        learning_rate=args.lr,
        weight_decay=args.weight_decay,
        warmup_ratio=args.warmup_ratio,
        eval_strategy="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        metric_for_best_model="f1_macro",
        greater_is_better=True,
        logging_dir=str(output_dir / "logs"),
        logging_steps=50,
        report_to=["none"],
        seed=args.seed,
        dataloader_num_workers=4,
        remove_unused_columns=False,
    )

    # Data collator: pad input_ids and attention_mask
    # Geneformer pad_token_id is 0 (<pad>)
    class GeneformerCollator:
        def __init__(self, pad_token_id=0):
            self.pad_token_id = pad_token_id

        def __call__(self, features):
            max_len = max(len(f["input_ids"]) for f in features)
            batch = {
                "input_ids": [],
                "attention_mask": [],
                "labels": [],
            }
            for f in features:
                ids = f["input_ids"]
                mask = f["attention_mask"]
                pad_len = max_len - len(ids)
                batch["input_ids"].append(ids + [self.pad_token_id] * pad_len)
                batch["attention_mask"].append(mask + [0] * pad_len)
                batch["labels"].append(f["labels"])

            return {
                k: torch.tensor(v) for k, v in batch.items()
            }

    collator = GeneformerCollator(pad_token_id=0)

    # Trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_ds,
        eval_dataset=val_ds,
        data_collator=collator,
        compute_metrics=compute_metrics,
        callbacks=[EarlyStoppingCallback(early_stopping_patience=2)],
    )

    print("Starting training...")
    trainer.train()

    # Evaluate on test set
    print("\nEvaluating on test set...")
    test_results = trainer.evaluate(test_ds)
    print(f"Test results: {test_results}")

    # Save model
    print(f"Saving model to {output_dir}")
    trainer.save_model(output_dir)

    # Save label mapping
    with open(output_dir / "label_mapping.json", "w") as f:
        json.dump({"label2id": label2id, "id2label": id2label}, f, indent=2)

    # Save test results
    with open(output_dir / "test_metrics.json", "w") as f:
        json.dump(test_results, f, indent=2)

    print("Done!")


if __name__ == "__main__":
    main()
