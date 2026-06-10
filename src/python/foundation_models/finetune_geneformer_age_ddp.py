"""
Fine-tune Geneformer V2-104M for age group classification with DDP (4 GPUs).

Same logic as single-GPU version but uses torchrun for DistributedDataParallel.
"""

import argparse
import json
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from datasets import load_from_disk
from sklearn.model_selection import train_test_split
from transformers import (
    BertForSequenceClassification,
    TrainingArguments,
    Trainer,
    EarlyStoppingCallback,
)


def freeze_layers(model, n_freeze):
    for param in model.bert.embeddings.parameters():
        param.requires_grad = False
    for i, layer in enumerate(model.bert.encoder.layer):
        if i < n_freeze:
            for param in layer.parameters():
                param.requires_grad = False


def compute_metrics(eval_pred):
    from sklearn.metrics import accuracy_score, f1_score
    logits, labels = eval_pred
    preds = np.argmax(logits, axis=1)
    return {
        "accuracy": accuracy_score(labels, preds),
        "f1_macro": f1_score(labels, preds, average="macro"),
        "f1_weighted": f1_score(labels, preds, average="weighted"),
    }


class GeneformerCollator:
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/tokenized/aida_v1_tokenized.dataset")
    parser.add_argument("--model-path", default="/home/scroll/zzhang/transcriptome-aging/data/geneformer/weights/geneformer-v2/Geneformer-V2-104M")
    parser.add_argument("--output-dir", default="/home/scroll/zzhang/transcriptome-aging/models/geneformer_age_classifier_ddp")
    parser.add_argument("--freeze-layers", type=int, default=6)
    parser.add_argument("--epochs", type=int, default=3)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--lr", type=float, default=5e-5)
    parser.add_argument("--weight-decay", type=float, default=0.01)
    parser.add_argument("--warmup-ratio", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # DDP setup
    local_rank = int(os.environ.get("LOCAL_RANK", -1))
    if local_rank == -1:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        torch.cuda.set_device(local_rank)
        device = torch.device("cuda", local_rank)
        torch.distributed.init_process_group(backend="nccl")

    if local_rank <= 0:
        print(f"Loading dataset...")
    ds = load_from_disk(args.dataset)

    label_list = sorted(set(ds["age_group"]))
    label2id = {l: i for i, l in enumerate(label_list)}
    id2label = {i: l for l, i in label2id.items()}

    def preprocess(example):
        example["labels"] = label2id[example["age_group"]]
        example["attention_mask"] = [1] * len(example["input_ids"])
        return example

    ds = ds.map(preprocess, num_proc=4)

    # Stratified split by donor
    donors = np.array(ds["donor_id"])
    labels = np.array(ds["labels"])
    donor_df = pd.DataFrame({"donor": donors, "label": labels})
    donor_labels = donor_df.groupby("donor")["label"].first().reset_index()

    train_donors, temp_donors = train_test_split(
        donor_labels["donor"].values, test_size=0.3,
        stratify=donor_labels["label"].values, random_state=args.seed,
    )
    val_donors, test_donors = train_test_split(
        temp_donors, test_size=0.5,
        stratify=donor_labels[donor_labels["donor"].isin(temp_donors)]["label"].values,
        random_state=args.seed,
    )

    train_mask = np.isin(donors, train_donors)
    val_mask = np.isin(donors, val_donors)
    test_mask = np.isin(donors, test_donors)

    train_ds = ds.select(np.where(train_mask)[0])
    val_ds = ds.select(np.where(val_mask)[0])
    test_ds = ds.select(np.where(test_mask)[0])

    if local_rank <= 0:
        print(f"Train: {len(train_ds)}, Val: {len(val_ds)}, Test: {len(test_ds)}")

    if local_rank <= 0:
        print("Loading model...")
    model = BertForSequenceClassification.from_pretrained(
        args.model_path, num_labels=len(label_list),
        id2label=id2label, label2id=label2id,
    )
    if args.freeze_layers > 0:
        freeze_layers(model, args.freeze_layers)

    if local_rank <= 0:
        trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
        total = sum(p.numel() for p in model.parameters())
        print(f"Trainable: {trainable/1e6:.1f}M / {total/1e6:.1f}M")

    output_dir = Path(args.output_dir)
    if local_rank <= 0:
        output_dir.mkdir(parents=True, exist_ok=True)

    training_args = TrainingArguments(
        output_dir=str(output_dir),
        num_train_epochs=args.epochs,
        per_device_train_batch_size=args.batch_size,
        per_device_eval_batch_size=args.batch_size * 2,
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
        dataloader_num_workers=2,
        remove_unused_columns=False,
        fp16=True,
        ddp_find_unused_parameters=False,
    )

    collator = GeneformerCollator(pad_token_id=0)

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_ds,
        eval_dataset=val_ds,
        data_collator=collator,
        compute_metrics=compute_metrics,
        callbacks=[EarlyStoppingCallback(early_stopping_patience=2)],
    )

    if local_rank <= 0:
        print("Starting training...")
    trainer.train()

    if local_rank <= 0:
        print("\nEvaluating on test set...")
        test_results = trainer.evaluate(test_ds)
        print(f"Test results: {test_results}")
        trainer.save_model(output_dir)
        with open(output_dir / "label_mapping.json", "w") as f:
            json.dump({"label2id": label2id, "id2label": id2label}, f, indent=2)
        with open(output_dir / "test_metrics.json", "w") as f:
            json.dump(test_results, f, indent=2)
        print("Done!")

    if local_rank != -1:
        torch.distributed.destroy_process_group()


if __name__ == "__main__":
    import os
    main()
