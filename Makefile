# HIV-DR-Watch developer shortcuts.
#
# Convention: targets that touch the host env prefix `run-`; targets
# that touch the compose stack prefix `up-` / `down-`. Bare targets
# are "do the obvious thing" for a fresh checkout.

SHELL := /usr/bin/env bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c

DATA_DIR      ?= data
LANL_DIR      ?= $(DATA_DIR)/raw/lanl
INTERIM_DIR   ?= $(DATA_DIR)/interim
GOLD_DIR      ?= $(DATA_DIR)/gold
MODELS_DIR    ?= $(DATA_DIR)/models

INGEST_OUT    := $(INTERIM_DIR)/east_africa_pol.parquet
CURATED_OUT   := $(INTERIM_DIR)/curated/aligned.parquet
ALIGNED_FASTA := $(INTERIM_DIR)/curated/aligned.fasta
CV_JSON       := $(GOLD_DIR)/cluster_validation.json
BENCH_DIR     := benchmarks/results

.PHONY: help
help:
	@echo "Targets:"
	@echo "  make setup            # copy .env.example → .env if missing"
	@echo "  make up               # docker compose up -d"
	@echo "  make down             # docker compose down"
	@echo "  make logs             # follow compose logs"
	@echo "  make ps               # list running services"
	@echo "  make images           # build MAFFT + tn93 + IQ-TREE + HIV-TRACE images"
	@echo "  make worker           # start a Prefect worker on the host"
	@echo "  make ingest           # LANL-only ingest → bronze parquet"
	@echo "                        # append TO_POSTGRES=1 to also load into DB"
	@echo "  make sierra           # sierra DR calls (cache-first)"
	@echo "  make curate           # align + trim + QC + dedup → silver"
	@echo "  make cluster          # tn93 + assignments + outbreak features"
	@echo "  make cluster-validate # our clusterer vs ground truth"
	@echo "                        # add WITH_HIVTRACE=1 for parity vs HIV-TRACE"
	@echo "  make phylo            # per-subtype IQ-TREE + TreeTime"
	@echo "  make phylo-smoke      # phylo, subtype C only (fast)"
	@echo "  make fig1             # pipeline architecture diagram"
	@echo "  make fig2             # resistance landscape heatmap"
	@echo "  make fig3             # cluster-validation + parity vs HIV-TRACE"
	@echo "  make fig4             # DTG case study (preregistered)"
	@echo "  make fig5             # scalability plot"
	@echo "  make paper1           # regenerate all 5 Paper 1 figures"
	@echo "  make benchmarks       # run every scale_* benchmark"
	@echo "  make benchmark-ingest / -align / -cluster"
	@echo "  make test             # pytest"

# ---------- setup ----------

.env:
	@cp .env.example .env
	@echo "wrote .env — edit if you need to remap ports"

.PHONY: setup
setup: .env
	uv sync

# ---------- compose stack ----------

.PHONY: up down logs ps
up: .env
	docker compose up -d
	@echo "Prefect: http://localhost:$${PREFECT_PORT:-4200}   MLflow: http://localhost:$${MLFLOW_PORT:-5000}   MinIO: http://localhost:$${MINIO_CONSOLE_PORT:-9001}"

down:
	docker compose down

logs:
	docker compose logs -f

ps:
	docker compose ps

# ---------- batch container images ----------

.PHONY: images
images:
	docker build -t hivdrwatch/mafft:latest docker/mafft
	docker build -t hivdrwatch/tn93:1.0.17 docker/hivtrace
	docker build -t hivdrwatch/hivtrace:0.1 docker/hivtrace_full
	docker build -t hivdrwatch/phylo:latest docker/iqtree

# ---------- Paper 1 figures ----------

.PHONY: fig1
fig1:
	python notebooks/fig01_pipeline_architecture.py

.PHONY: fig2
fig2:
	python notebooks/fig02_resistance_landscape.py

.PHONY: fig3
fig3:
	@test -f $(CV_JSON) || { echo "$(CV_JSON) missing — run 'make cluster-validate WITH_HIVTRACE=1' first"; exit 1; }
	python notebooks/fig03_cluster_validation.py

.PHONY: fig4
fig4:
	python notebooks/fig04_dtg_case_study.py

.PHONY: fig5
fig5:
	@ls $(BENCH_DIR)/scale_*.json >/dev/null 2>&1 || { echo "$(BENCH_DIR)/scale_*.json missing — run 'make benchmarks' first"; exit 1; }
	python notebooks/fig05_scalability.py

.PHONY: paper1
paper1: fig1 fig2 fig3 fig4 fig5
	@echo "Paper 1 figures regenerated under paper1_tools/figures/"

# ---------- benchmarks (feed Fig 5, checked into git for CI diff) ----------

.PHONY: benchmarks benchmark-ingest benchmark-align benchmark-cluster
benchmarks: benchmark-ingest benchmark-align benchmark-cluster

benchmark-ingest:
	python -m benchmarks.scale_ingest \
	    --sizes $${SIZES:-100 500 1000 5000} \
	    --out $(BENCH_DIR)/scale_ingest.json \
	    --seed $${SEED:-42}

benchmark-align:
	python -m benchmarks.scale_align \
	    --sizes $${SIZES:-100 200 500 1000 2000} \
	    --out $(BENCH_DIR)/scale_align.json \
	    --seed $${SEED:-42}

benchmark-cluster:
	python -m benchmarks.scale_cluster \
	    --sizes $${SIZES:-100 200 500} \
	    --out $(BENCH_DIR)/scale_cluster.json \
	    --seed $${SEED:-42}

# ---------- cluster validation (Fig 3) ----------

.PHONY: cluster-validate
cluster-validate:
	python -m hivdrwatch.pipeline.cluster_validate \
	    --out $(CV_JSON) \
	    --work-dir $(INTERIM_DIR)/cluster_validate \
	    --n-clusters $${N_CLUSTERS:-20} \
	    --cluster-size-lo $${CLUSTER_SIZE_LO:-5} \
	    --cluster-size-hi $${CLUSTER_SIZE_HI:-15} \
	    --n-singletons $${N_SINGLETONS:-50} \
	    --seed $${SEED:-42} \
	    $${WITH_HIVTRACE:+--with-hivtrace}

.PHONY: phylo
phylo:
	python -m hivdrwatch.pipeline.phylo \
	    --curated $(CURATED_OUT) \
	    --trees-root $(INTERIM_DIR)/trees \
	    --manifest data/gold/trees/manifest.parquet \
	    --seed $${SEED:-42} \
	    --bootstrap $${BOOTSTRAP:-0} \
	    --min-coverage $${MIN_COVERAGE:-0.5} \
	    --min-taxa $${MIN_TAXA:-20} \
	    $${NO_FAST:+--no-fast} \
	    $${SUBTYPES:+--subtypes $$SUBTYPES}

# Fast smoke of the whole phylo path — builds only the smallest subtype
# (C, ~20 seqs) to prove the container + wrapper + flow all work.
.PHONY: phylo-smoke
phylo-smoke:
	SUBTYPES=C $(MAKE) phylo

# ---------- host worker ----------

.PHONY: worker
worker:
	prefect worker start --pool default-agent-pool

.PHONY: restart-worker
restart-worker:
	@# Kill any running prefect worker; ok if none found.
	@pkill -f "prefect worker start" 2>/dev/null || true
	@sleep 1
	@echo "old workers killed; starting fresh worker..."
	@prefect worker start --pool default-agent-pool &
	@sleep 2
	@echo "worker restarted with current code"

# ---------- pipeline stages ----------

.PHONY: ingest
ingest:
	python -m hivdrwatch.ingest --skip-genbank --lanl-dir $(LANL_DIR) --out $(INGEST_OUT) $${TO_POSTGRES:+--to-postgres}

.PHONY: ingest-genbank
ingest-genbank:
	@[ -n "$${NCBI_EMAIL}" ] || (echo "set NCBI_EMAIL in .env before ingest-genbank" && exit 1)
	python -m hivdrwatch.pipeline.ingest \
	    --email $${NCBI_EMAIL} \
	    --countries KE UG TZ RW BI ET \
	    --lanl-dir $(LANL_DIR) \
	    --genbank-dir $(DATA_DIR)/raw/genbank \
	    --out $(INGEST_OUT) \
	    $${TO_POSTGRES:+--to-postgres}

.PHONY: sierra
sierra:
	python -m hivdrwatch.pipeline.dr_ingest --input $(INGEST_OUT) --cache-dir $(DATA_DIR)/raw/hivdb --out-dir $(INTERIM_DIR) $${TO_POSTGRES:+--to-postgres}

.PHONY: curate
curate:
	python -m hivdrwatch.pipeline.curate --input $(INGEST_OUT) --lanl-dir $(LANL_DIR) --work-dir $(INTERIM_DIR)/curate --out $(CURATED_OUT) --aligner $${ALIGNER:-mappy} --chunk-size $${CHUNK_SIZE:-1000} --chunk-concurrency $${CHUNK_CONCURRENCY:-2}

.PHONY: cluster
cluster:
	python -m hivdrwatch.pipeline.cluster --aligned-fasta $(ALIGNED_FASTA) --curated $(CURATED_OUT) --work-dir $(INTERIM_DIR)/cluster --gold-dir $(GOLD_DIR) $${TO_POSTGRES:+--to-postgres}

.PHONY: pipeline
pipeline: ingest sierra curate cluster
	@echo "done — curated + DR + clusters under $(INTERIM_DIR) and $(GOLD_DIR)"

# ---------- storage (Batch E) ----------

.PHONY: storage-migrate
storage-migrate:
	alembic -c src/hivdrwatch/storage/alembic.ini upgrade head

.PHONY: storage-load-sequences
storage-load-sequences:
	python -c "from pathlib import Path; from hivdrwatch.storage.session import get_session; from hivdrwatch.storage.loader import load_sequences, ensure_schema; \
	sess_ctx = get_session(); sess = sess_ctx.__enter__(); ensure_schema(sess); counts = load_sequences(Path('$(INGEST_OUT)'), session=sess); print(counts); sess_ctx.__exit__(None, None, None)"

# ---------- api (Batch F) ----------

.PHONY: api
api:
	docker compose up -d api
	@echo "API on http://localhost:$${API_PORT:-8080}   docs at /docs"

.PHONY: api-dev
api-dev:
	uvicorn hivdrwatch.api.main:app --host 0.0.0.0 --port $${API_PORT:-8080} --reload

# ---------- foundation_models / Paper 3 baselines ----------

# Fit the zero-GPU baselines (k-mer LR + one-hot RF) against Paper 1
# outputs and write per-drug × per-subtype metrics. GPU-required
# ESM-2 / Evo fine-tuning lands separately.
.PHONY: fm-baselines
fm-baselines:
	python -c "\
import polars as pl; \
from pathlib import Path; \
from foundation_models import build_fm_dataset; \
from foundation_models.baselines import KmerLogisticRegression; \
from foundation_models.evaluate import per_subtype_metrics; \
dr = pl.read_parquet('$(INTERIM_DIR)/dr_results.parquet'); \
cur = pl.read_parquet('$(CURATED_OUT)'); \
ds = build_fm_dataset(dr, cur, strategy='$${STRATEGY:-subtype_holdout}', seed=42); \
rows = []; \
for drug, test_df in ds.test.partition_by('drug', as_dict=True).items(): \
    train_df = ds.train.filter(pl.col('drug') == drug); \
    if train_df.is_empty(): \
        continue; \
    m = KmerLogisticRegression(k=3, seed=42).fit(train_df['sequence'].to_list(), train_df['label_binary'].to_list()); \
    probs = m.predict_proba(test_df['sequence'].to_list()); \
    rows.append(per_subtype_metrics(test_df.with_columns(pl.Series('p', probs)), prob_col='p')); \
out = Path('$(GOLD_DIR)/paper3_baseline_metrics.parquet'); \
out.parent.mkdir(parents=True, exist_ok=True); \
(pl.concat(rows) if rows else pl.DataFrame()).write_parquet(out); \
print(f'wrote {out}')"

# ---------- causal / Paper 2 ----------

.PHONY: causal
causal:
	python -m causal.flows.causal \
	    --dr-parquet $(INTERIM_DIR)/dr_results.parquet \
	    --curated-parquet $(CURATED_OUT) \
	    --out-dir $(GOLD_DIR) \
	    --bootstrap-samples $${BOOTSTRAP:-200} \
	    $${PER_SUBTYPE:+--per-subtype} \
	    $${PER_DRUG:+--per-drug}

# ---------- paper 2 (CausalPhylo) manuscript ----------

.PHONY: paper2 paper2-bench paper2-fig1 paper2-render
paper2: paper2-bench paper2-fig1 paper2-render

# Synthetic-ground-truth benchmark that generates the numbers cited
# in paper2_causal/manuscript.qmd §Results.
paper2-bench:
	python -m benchmarks.paper2_synth_recovery \
	    --sizes 500 2000 8000 \
	    --seed 42 \
	    --out benchmarks/results/paper2_synth_recovery.json

# Render the DR domain DAG figure from causal/dag.py — no external
# data needed.
paper2-fig1:
	python paper2_causal/notebooks/paper2_fig01_dag.py

# Render the Quarto manuscript to PDF + HTML.
paper2-render:
	cd paper2_causal && quarto render manuscript.qmd

# ---------- paper 3 (foundation-model benchmark) manuscript ----------

.PHONY: paper3 paper3-bench paper3-render
paper3: paper3-bench paper3-render

paper3-bench:
	python -m benchmarks.paper3_baselines_synth \
	    --sizes 2000 8000 20000 \
	    --seed 42 \
	    --out benchmarks/results/paper3_baselines_synth.json

paper3-render:
	cd paper3_foundation && quarto render manuscript.qmd

# ---------- EDA + supplement notebooks (Marimo) ----------

.PHONY: eda dtg-supp
# Comprehensive v0.2.0 exploratory data-analysis report — A/B/C/D/E/F/H/I.
# Renders to a self-contained HTML embeddable in the SAWA dashboard.
eda:
	marimo export html notebooks/eda_v0.2.0.py \
	    -o paper1_tools/figures/eda_v0.2.0.html

# Fig 4 supplement — DTG-era INSTI resistance deep-dive (§G).
# Per-position INSTI mutations, per-drug level distributions,
# preregistered sensitivity analyses, BH-FDR correction diagnostic.
dtg-supp:
	marimo export html notebooks/dtg_deep_dive.py \
	    -o paper1_tools/figures/dtg_deep_dive.html

# ---------- DVC (data provenance) ----------

.PHONY: dvc-init dvc-track dvc-status
# Idempotent DVC bootstrap: init if needed, then track every gold-tier
# artifact. Safe to re-run after new data lands; only adds pointers.
dvc-init dvc-track:
	bash scripts/setup_dvc.sh

# Quick check: what's out-of-sync between tracked pointers and disk?
dvc-status:
	dvc status

# ---------- Foundation-model training (Paper 3) ----------

.PHONY: fm-prepare fm-esm-smoke fm-esm-train fm-esm-all
# Build train/val/test parquets per drug class from the curated + Sierra
# outputs. Writes to data/gold/fm/<class>/<strategy>/ plus a manifest.
fm-prepare:
	python -m foundation_models.prepare_training_data \
	    --curated data/interim/curated/aligned.parquet \
	    --dr-results data/interim/dr_results.parquet \
	    --out-dir data/gold/fm \
	    --drug-classes INSTI NRTI NNRTI PI \
	    --strategy stratified \
	    --seed 42

# ---- Shared FM training hyperparameters ----
# EPOCHS, BATCH_SIZE, GRAD_ACCUM live here and ONLY here. Every ESM-2
# training target below must reference them so stratified and
# subtype_holdout runs stay experimentally comparable — the whole
# point of the two splits is to hold everything else constant and
# vary only the train/test partition. If you change EPOCHS to compare
# runs at different training lengths, override on the command line
# for both invocations, not by patching one target.
#
#   make fm-esm-train EPOCHS=20 CLASS=INSTI
#   make fm-esm-holdout EPOCHS=20 CLASS=INSTI   # SAME EPOCHS
EPOCHS      ?= 10
BATCH_SIZE  ?= 16
GRAD_ACCUM  ?= 2
CLASS       ?= INSTI
STRATEGY    ?= stratified
SEED        ?= 42

# 1-epoch, 500-row smoke test on INSTI. Verifies GPU + tokenizer +
# training loop plumbing before you commit compute to a real run.
# Expected wall: <60s on any GPU with 4+GB VRAM. No MLflow.
fm-esm-smoke:
	python -m foundation_models.esm \
	    --drug-class INSTI --strategy stratified --smoke --verbose

# Fine-tune on the stratified split. Uses shared EPOCHS/BATCH_SIZE so
# fm-esm-holdout is a proper controlled comparison.
fm-esm-train:
	python -m foundation_models.esm \
	    --drug-class $(CLASS) --strategy stratified \
	    --batch-size $(BATCH_SIZE) --grad-accum-steps $(GRAD_ACCUM) \
	    --max-epochs $(EPOCHS) --seed $(SEED) \
	    $(if $(MLFLOW_URI),--mlflow-uri $(MLFLOW_URI),)

# All four drug classes back-to-back on the stratified split.
fm-esm-all: fm-prepare
	for cls in INSTI NRTI NNRTI PI; do \
	  python -m foundation_models.esm \
	      --drug-class $$cls --strategy stratified \
	      --batch-size $(BATCH_SIZE) --grad-accum-steps $(GRAD_ACCUM) \
	      --max-epochs $(EPOCHS) --seed $(SEED) \
	      $(if $(MLFLOW_URI),--mlflow-uri $(MLFLOW_URI),); \
	done

# Head-to-head vs kmer-LR + one-hot-RF baselines on the SAME test split
# the ESM-2 checkpoint was scored against. This is the essential Paper 3
# comparison — ESM-2 test AUROC is only meaningful next to a baseline
# on identical data. Emits data/gold/fm/comparison_<class>_<strategy>.{parquet,json}.
#
# Vars: CLASS=INSTI (default), STRATEGY=stratified (default).
# Example: make fm-compare CLASS=NRTI STRATEGY=subtype_holdout
fm-compare:
	python -m foundation_models.compare_baselines \
	    --drug-class $(CLASS) --strategy $(STRATEGY) --verbose

# Subtype-holdout variant of the ESM-2 fine-tune — Paper 3's thesis-
# testing split. Train on {A1, C, D} (dominant East African), test on
# {B, unknown, recombinant, minor subtypes}. Uses the SAME shared
# EPOCHS/BATCH_SIZE/SEED as fm-esm-train so the stratified vs holdout
# comparison isolates only the train/test partition.
fm-esm-holdout:
	python -m foundation_models.prepare_training_data \
	    --curated data/interim/curated/aligned.parquet \
	    --dr-results data/interim/dr_results.parquet \
	    --out-dir data/gold/fm \
	    --drug-classes INSTI NRTI NNRTI PI \
	    --strategy subtype_holdout \
	    --holdout-train-subtypes A1 C D \
	    --seed $(SEED)
	python -m foundation_models.esm \
	    --drug-class $(CLASS) --strategy subtype_holdout \
	    --manifest data/gold/fm/manifest_subtype_holdout.json \
	    --batch-size $(BATCH_SIZE) --grad-accum-steps $(GRAD_ACCUM) \
	    --max-epochs $(EPOCHS) --seed $(SEED) \
	    $(if $(MLFLOW_URI),--mlflow-uri $(MLFLOW_URI),)

# ---------- Problems log (methodological record) ----------

.PHONY: log-problem
# Interactive prompt to append a new PL-### entry to docs/PROBLEMS_LOG.md.
# Auto-assigns the next ID, records today's date + current git SHA,
# refuses to run if PROBLEMS_LOG.md has uncommitted edits (pass --force
# to override). Every future methodological issue gets logged this way.
log-problem:
	python scripts/log_problem.py

# ---------- tests ----------

.PHONY: test
test:
	pytest tests/ -q
