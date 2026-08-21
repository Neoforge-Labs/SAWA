"""ESM-2 (protein LM) fine-tuning scaffold for HIV-1 drug-resistance.

Default backbone: ``facebook/esm2_t6_8M_UR50D`` — 6 transformer layers,
320-dim hidden, ~8M parameters. Chosen because:

* Fits comfortably in <2 GB VRAM at batch_size=32 for pol sequences
  (~830 AA + special tokens ≈ 830 seq positions), so it runs on a
  mid-tier consumer GPU.
* Small enough to fully fine-tune (no LoRA required); large enough to
  meaningfully improve over kmer-LR / one-hot-RF baselines if there's
  signal beyond n-grams.
* Same tokenizer family as the larger ESM-2 checkpoints, so a
  scale-up experiment (150M / 650M) is a one-line change.

**All heavy imports (torch, transformers, mlflow) are deferred inside
functions.** The rest of the codebase must remain importable without
these deps — install with ``uv sync --extra fm``.

Usage (from a script or the CLI at __main__):

    from foundation_models.esm import ESM2Config, train_from_manifest

    cfg = ESM2Config(
        drug_class="INSTI",
        strategy="stratified",
        max_epochs=3,
        batch_size=16,
    )
    train_from_manifest(
        manifest="data/gold/fm/manifest_stratified.json",
        cfg=cfg,
        ckpt_dir="data/gold/fm/checkpoints",
    )

Or from the terminal (see ``__main__`` at the bottom):

    python -m foundation_models.esm \\
        --drug-class INSTI --strategy stratified \\
        --smoke  # 1 epoch on 500 rows, no MLflow
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

DEFAULT_MODEL = "facebook/esm2_t6_8M_UR50D"
DEFAULT_MAX_LENGTH = 1024  # ESM-2 max positional embedding
POL_MAX_AA = 850  # PR + RT + IN typical


@dataclass
class ESM2Config:
    """Training configuration — everything reproducible from these fields."""

    drug_class: str
    strategy: str = "stratified"
    model_name: str = DEFAULT_MODEL
    max_length: int = POL_MAX_AA
    batch_size: int = 16
    grad_accum_steps: int = 2
    learning_rate: float = 2e-5
    weight_decay: float = 0.01
    max_epochs: int = 3
    warmup_steps: int = 50
    seed: int = 42
    fp16: bool = True
    freeze_backbone: bool = False
    early_stop_patience: int = 2
    positive_class_weight: float | None = None  # auto-computed if None
    log_every_n_steps: int = 20
    device: str = "auto"  # "auto" | "cuda" | "cpu"
    smoke: bool = False  # if True, cap train to 500 rows, 1 epoch

    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class TrainResult:
    """Everything a downstream evaluator needs to compare runs."""

    drug_class: str
    strategy: str
    n_train: int
    n_val: int
    n_test: int
    best_val_auroc: float
    test_auroc: float
    test_auprc: float
    wall_seconds: float
    ckpt_path: str
    config: dict


# ------------------------------------------------------------------
# Torch / HuggingFace bits (deferred imports throughout)
# ------------------------------------------------------------------


def _resolve_device(pref: str) -> str:
    """Pick the best available torch device.

    On "auto":
      * CUDA if a discrete NVIDIA GPU is visible (fastest).
      * MPS if we're on Apple Silicon with a recent torch (M1/M2/M3
        Macs — the "Metal Performance Shaders" backend; ~5-10x faster
        than CPU on ESM-2 t6 in practice).
      * CPU fallback with a clear diagnostic — most common cause is a
        CPU-only torch wheel; the fix is a CUDA-specific or Mac wheel.
    """
    import torch

    if pref != "auto":
        return pref
    if torch.cuda.is_available():
        return "cuda"
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "mps"
    logger.warning(
        "no GPU visible to torch — falling back to CPU. Diagnose with:\n"
        "  python -c 'import torch; "
        'print(dict(cuda=torch.cuda.is_available(), '
        'mps=getattr(torch.backends, \\"mps\\", None) and torch.backends.mps.is_available()))'
        "'\n"
        "  If both False on a Mac, reinstall torch WITHOUT --index-url overrides.\n"
        "  On Linux/Windows with NVIDIA GPU: install the matching cu12x wheel."
    )
    return "cpu"


def _load_splits(manifest_path: Path, drug_class: str, strategy: str):
    """Load train/val/test parquets for a drug class from the FM manifest."""
    import polars as pl

    manifest = json.loads(manifest_path.read_text())
    runs = manifest.get("runs", [])
    match = next(
        (
            r
            for r in runs
            if r["drug_class"] == drug_class and r["strategy"] == strategy
        ),
        None,
    )
    if match is None:
        raise ValueError(
            f"manifest {manifest_path} has no entry for "
            f"drug_class={drug_class!r} strategy={strategy!r}. "
            f"Available: {[(r['drug_class'], r['strategy']) for r in runs]}"
        )
    train = pl.read_parquet(match["train_parquet"])
    val = pl.read_parquet(match["val_parquet"])
    test = pl.read_parquet(match["test_parquet"])
    return train, val, test, match


class _SequenceLabelDataset:
    """Torch Dataset over (sequence_aa, label_binary) rows."""

    def __init__(self, df, tokenizer, max_length: int) -> None:
        from torch.utils.data import Dataset  # noqa: F401

        # Store as parallel python lists — Polars → list conversion is
        # a one-time cost, and downstream iteration is fast.
        self.sequences = df["sequence_aa"].to_list()
        self.labels = df["label_binary"].to_list()
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, idx: int) -> dict:
        enc = self.tokenizer(
            self.sequences[idx],
            truncation=True,
            max_length=self.max_length,
            padding=False,
            return_tensors=None,
        )
        enc["labels"] = int(self.labels[idx])
        return enc


def _make_collator(tokenizer):
    from transformers import DataCollatorWithPadding

    return DataCollatorWithPadding(tokenizer=tokenizer)


def _label_balanced_subsample(df, n: int, seed: int):
    """Take ``n`` rows preserving the source frame's positive rate.

    Falls back to a plain shuffle-and-head if there are no positives
    to preserve. Deterministic given seed.
    """
    import polars as pl

    n_pos = int(df["label_binary"].sum())
    if n_pos == 0 or df.height <= n:
        return df.sample(n=min(n, df.height), seed=seed) if df.height else df
    pos = df.filter(pl.col("label_binary") == 1)
    neg = df.filter(pl.col("label_binary") == 0)
    pos_frac = min(1.0, pos.height / max(1, df.height))
    n_pos_take = max(1, int(n * pos_frac))
    n_neg_take = max(0, n - n_pos_take)
    return pl.concat(
        [
            pos.sample(n=min(n_pos_take, pos.height), seed=seed),
            neg.sample(n=min(n_neg_take, neg.height), seed=seed),
        ]
    ).sample(fraction=1.0, seed=seed)  # shuffle so batches aren't sorted


def _compute_pos_weight(train_df) -> float:
    import torch

    n_pos = int(train_df["label_binary"].sum())
    n_neg = train_df.height - n_pos
    if n_pos == 0:
        logger.warning("no positive labels in train — pos_weight defaulting to 1.0")
        return 1.0
    return float(n_neg / n_pos)


def _metrics(labels, probs) -> tuple[float, float]:
    """Return (AUROC, AUPRC). Both guarded for degenerate label sets."""
    import numpy as np
    from sklearn.metrics import average_precision_score, roc_auc_score

    labels = np.asarray(labels)
    probs = np.asarray(probs)
    if labels.sum() == 0 or labels.sum() == labels.size:
        return (float("nan"), float("nan"))
    return float(roc_auc_score(labels, probs)), float(
        average_precision_score(labels, probs)
    )


def _evaluate(model, loader, device) -> tuple[float, float, list[int], list[float]]:
    import torch

    model.eval()
    all_labels: list[int] = []
    all_probs: list[float] = []
    with torch.no_grad():
        for batch in loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            labels = batch.pop("labels")
            logits = model(**batch).logits
            probs = torch.softmax(logits, dim=-1)[:, 1]
            all_labels.extend(labels.cpu().tolist())
            all_probs.extend(probs.float().cpu().tolist())
    auroc, auprc = _metrics(all_labels, all_probs)
    return auroc, auprc, all_labels, all_probs


def train_from_manifest(
    manifest: str | Path,
    cfg: ESM2Config,
    ckpt_dir: str | Path,
    *,
    mlflow_tracking_uri: str | None = None,
) -> TrainResult:
    """Fine-tune ESM-2 on the requested drug-class × strategy split.

    Reads train/val/test Parquets from the manifest emitted by
    ``foundation_models.prepare_training_data``. Writes the best-val
    checkpoint under
    ``<ckpt_dir>/<drug_class>/<strategy>/best.pt`` plus a
    ``run_metadata.json`` beside it.
    """
    # Deferred imports keep the module importable without torch.
    import numpy as np
    import torch
    from torch.optim import AdamW
    from torch.utils.data import DataLoader
    from transformers import (
        AutoModelForSequenceClassification,
        AutoTokenizer,
        get_linear_schedule_with_warmup,
    )

    manifest_path = Path(manifest)
    train_df, val_df, test_df, run_meta = _load_splits(
        manifest_path, cfg.drug_class, cfg.strategy
    )
    if cfg.smoke:
        # Preserve positive:negative balance — plain head() on Polars
        # returns natural sort order which is all-negatives-first for
        # ~2% positive-rate INSTI, wrecking val/test AUROC.
        train_df = _label_balanced_subsample(train_df, n=500, seed=cfg.seed)
        val_df = _label_balanced_subsample(val_df, n=200, seed=cfg.seed)
        test_df = _label_balanced_subsample(test_df, n=200, seed=cfg.seed)
        cfg.max_epochs = 1

    device = _resolve_device(cfg.device)
    logger.info(
        "device=%s | train=%d val=%d test=%d | model=%s | fp16=%s",
        device,
        train_df.height,
        val_df.height,
        test_df.height,
        cfg.model_name,
        cfg.fp16,
    )

    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)

    # MPS fp16 autocast silently triggers CPU fallback for many ops in
    # PyTorch 2.x on Apple Silicon — measured ~150s/step on M-series
    # for pol-length sequences vs ~2s in fp32. Override the config
    # default rather than surprise the operator.
    if device == "mps" and cfg.fp16:
        logger.warning(
            "MPS + fp16 is slow in PyTorch 2.x (silent CPU fallback for "
            "several ops). Overriding fp16=False for this run — pass "
            "--no-fp16 explicitly to silence this warning, or override "
            "with cfg.fp16=True after benchmarking your PyTorch version."
        )
        cfg.fp16 = False

    tokenizer = AutoTokenizer.from_pretrained(cfg.model_name)
    # Request scaled-dot-product-attention when available. On MPS this
    # hits Metal-optimized kernels; on CPU it falls back cleanly.
    try:
        model = AutoModelForSequenceClassification.from_pretrained(
            cfg.model_name, num_labels=2, attn_implementation="sdpa"
        ).to(device)
    except (ValueError, TypeError):
        # Older transformers or ESM versions don't accept the kwarg —
        # fall back to eager attention without failing the run.
        model = AutoModelForSequenceClassification.from_pretrained(
            cfg.model_name, num_labels=2
        ).to(device)

    if cfg.freeze_backbone:
        for name, p in model.named_parameters():
            if not name.startswith("classifier"):
                p.requires_grad = False
        logger.info("froze backbone — only the classifier head trains")

    train_ds = _SequenceLabelDataset(train_df, tokenizer, cfg.max_length)
    val_ds = _SequenceLabelDataset(val_df, tokenizer, cfg.max_length)
    test_ds = _SequenceLabelDataset(test_df, tokenizer, cfg.max_length)

    collator = _make_collator(tokenizer)
    train_loader = DataLoader(
        train_ds,
        batch_size=cfg.batch_size,
        shuffle=True,
        collate_fn=collator,
        num_workers=0,
    )
    val_loader = DataLoader(
        val_ds, batch_size=cfg.batch_size, shuffle=False, collate_fn=collator
    )
    test_loader = DataLoader(
        test_ds, batch_size=cfg.batch_size, shuffle=False, collate_fn=collator
    )

    pos_weight = cfg.positive_class_weight or _compute_pos_weight(train_df)
    class_weights = torch.tensor([1.0, pos_weight], dtype=torch.float, device=device)
    loss_fn = torch.nn.CrossEntropyLoss(weight=class_weights)

    optimizer = AdamW(
        [p for p in model.parameters() if p.requires_grad],
        lr=cfg.learning_rate,
        weight_decay=cfg.weight_decay,
    )
    total_steps = max(1, cfg.max_epochs * len(train_loader) // cfg.grad_accum_steps)
    scheduler = get_linear_schedule_with_warmup(
        optimizer, num_warmup_steps=cfg.warmup_steps, num_training_steps=total_steps
    )
    # GradScaler + fp16 autocast is CUDA-only. MPS supports autocast
    # without a scaler (float16 on-device); CPU stays fp32.
    scaler = torch.cuda.amp.GradScaler() if (cfg.fp16 and device == "cuda") else None
    if cfg.fp16 and device == "mps":
        logger.info("MPS: using autocast(float16) without GradScaler")
    if cfg.fp16 and device == "cpu":
        logger.info("CPU: ignoring fp16 request (float16 CPU is slower than fp32)")

    # MLflow is optional — if unavailable, log to stderr only.
    mlf = None
    if mlflow_tracking_uri:
        try:
            import mlflow as _mlf

            _mlf.set_tracking_uri(mlflow_tracking_uri)
            _mlf.set_experiment(f"fm_esm2_{cfg.drug_class.lower()}")
            mlf = _mlf.start_run(run_name=f"{cfg.drug_class}_{cfg.strategy}")
            _mlf.log_params(asdict(cfg))
        except ImportError:
            logger.warning("mlflow not installed; skipping tracking")

    best_val_auroc = -1.0
    best_state: dict | None = None
    patience = cfg.early_stop_patience
    stalled = 0

    t0 = time.time()
    for epoch in range(cfg.max_epochs):
        model.train()
        optimizer.zero_grad()
        for step, batch in enumerate(train_loader):
            batch = {k: v.to(device) for k, v in batch.items()}
            labels = batch.pop("labels")

            if scaler is not None:
                with torch.cuda.amp.autocast(dtype=torch.float16):
                    logits = model(**batch).logits
                    loss = loss_fn(logits, labels) / cfg.grad_accum_steps
                scaler.scale(loss).backward()
                if (step + 1) % cfg.grad_accum_steps == 0:
                    scaler.step(optimizer)
                    scaler.update()
                    scheduler.step()
                    optimizer.zero_grad()
            elif cfg.fp16 and device == "mps":
                # MPS autocast — no scaler; MPS float16 is on-device.
                with torch.autocast(device_type="mps", dtype=torch.float16):
                    logits = model(**batch).logits
                    loss = loss_fn(logits, labels) / cfg.grad_accum_steps
                loss.backward()
                if (step + 1) % cfg.grad_accum_steps == 0:
                    optimizer.step()
                    scheduler.step()
                    optimizer.zero_grad()
            else:
                logits = model(**batch).logits
                loss = loss_fn(logits, labels) / cfg.grad_accum_steps
                loss.backward()
                if (step + 1) % cfg.grad_accum_steps == 0:
                    optimizer.step()
                    scheduler.step()
                    optimizer.zero_grad()

            if step % cfg.log_every_n_steps == 0:
                logger.info(
                    "epoch=%d step=%d loss=%.4f lr=%.2e",
                    epoch,
                    step,
                    float(loss.detach()) * cfg.grad_accum_steps,
                    scheduler.get_last_lr()[0],
                )

        val_auroc, val_auprc, _, _ = _evaluate(model, val_loader, device)
        logger.info(
            "epoch=%d done | val AUROC=%.4f val AUPRC=%.4f",
            epoch,
            val_auroc,
            val_auprc,
        )
        if mlf is not None:
            import mlflow

            mlflow.log_metrics(
                {"val_auroc": val_auroc, "val_auprc": val_auprc}, step=epoch
            )

        if val_auroc > best_val_auroc:
            best_val_auroc = val_auroc
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            stalled = 0
        else:
            stalled += 1
            if stalled >= patience:
                logger.info("early stopping at epoch %d (patience=%d)", epoch, patience)
                break

    # Load best, evaluate on test.
    if best_state is not None:
        model.load_state_dict(best_state)
    test_auroc, test_auprc, _, _ = _evaluate(model, test_loader, device)
    wall = time.time() - t0
    logger.info(
        "test AUROC=%.4f AUPRC=%.4f | wall=%.1fs", test_auroc, test_auprc, wall
    )

    ckpt_dir = Path(ckpt_dir) / cfg.drug_class.lower() / cfg.strategy
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    ckpt_path = ckpt_dir / "best.pt"
    torch.save(
        {
            "model_state_dict": best_state or model.state_dict(),
            "config": asdict(cfg),
            "best_val_auroc": best_val_auroc,
        },
        ckpt_path,
    )
    result = TrainResult(
        drug_class=cfg.drug_class,
        strategy=cfg.strategy,
        n_train=train_df.height,
        n_val=val_df.height,
        n_test=test_df.height,
        best_val_auroc=best_val_auroc,
        test_auroc=test_auroc,
        test_auprc=test_auprc,
        wall_seconds=wall,
        ckpt_path=str(ckpt_path),
        config=asdict(cfg),
    )
    (ckpt_dir / "run_metadata.json").write_text(
        json.dumps(asdict(result), indent=2), encoding="utf-8"
    )

    if mlf is not None:
        import mlflow

        mlflow.log_metrics(
            {"test_auroc": test_auroc, "test_auprc": test_auprc, "wall_s": wall}
        )
        mlflow.end_run()

    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="foundation_models.esm")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/gold/fm/manifest_stratified.json"),
    )
    parser.add_argument("--drug-class", required=True)
    parser.add_argument(
        "--strategy", choices=("stratified", "subtype_holdout"), default="stratified"
    )
    parser.add_argument("--model-name", default=DEFAULT_MODEL)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--grad-accum-steps", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=2e-5)
    parser.add_argument("--max-epochs", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--freeze-backbone", action="store_true")
    parser.add_argument("--no-fp16", action="store_true")
    parser.add_argument(
        "--smoke", action="store_true", help="1 epoch on 500 rows, no MLflow"
    )
    parser.add_argument("--device", default="auto")
    parser.add_argument(
        "--ckpt-dir", type=Path, default=Path("data/gold/fm/checkpoints")
    )
    parser.add_argument("--mlflow-uri", default=None)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cfg = ESM2Config(
        drug_class=args.drug_class,
        strategy=args.strategy,
        model_name=args.model_name,
        batch_size=args.batch_size,
        grad_accum_steps=args.grad_accum_steps,
        learning_rate=args.learning_rate,
        max_epochs=args.max_epochs,
        seed=args.seed,
        freeze_backbone=args.freeze_backbone,
        fp16=not args.no_fp16,
        smoke=args.smoke,
        device=args.device,
    )
    train_from_manifest(
        args.manifest,
        cfg=cfg,
        ckpt_dir=args.ckpt_dir,
        mlflow_tracking_uri=args.mlflow_uri if not args.smoke else None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
