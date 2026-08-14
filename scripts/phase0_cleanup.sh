#!/usr/bin/env bash
# Phase 0 cleanup — align the repo to MODULES.md contract.
#
# Executes the four decisions from the 2026-08-13 architect audit:
#   1. Delete out-of-scope modules (models/, features/, naomi/) + their flows/tests
#   2. Move causal/ from src/hivdrwatch/ to repo root
#   3. Split flows/ into pipeline/ (Prefect live) — keep in place, rename
#      (workflows/ for future Nextflow paper analyses is created empty)
#   4. Delete dashboard notebook
#   5. Delete data/gold/ + data/models/ + naomi labels
#
# Uses git mv where possible to preserve history. Run from repo root.
# Idempotent-ish: skips things that don't exist.
#
# Usage:
#   bash scripts/phase0_cleanup.sh          # dry-run: prints what it would do
#   bash scripts/phase0_cleanup.sh --apply  # actually do it

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1

run() {
    if [[ $APPLY -eq 1 ]]; then
        echo "  $*"
        eval "$@"
    else
        echo "  [dry-run] $*"
    fi
}

# ------------------------------------------------------------------
# Batch 1: delete out-of-scope modules + tests
# ------------------------------------------------------------------

echo ""
echo "== Batch 1: delete out-of-scope modules =="
for d in \
    src/hivdrwatch/models \
    src/hivdrwatch/features \
    src/hivdrwatch/naomi \
    tests/models \
    tests/features \
    tests/naomi
do
    if [[ -e "$d" ]]; then
        run "git rm -rf '$d' || rm -rf '$d'"
    fi
done

# Individual flow files that supported deleted modules
for f in \
    src/hivdrwatch/flows/dr_features.py \
    src/hivdrwatch/flows/train_dr.py \
    src/hivdrwatch/flows/train_outbreak.py \
    src/hivdrwatch/flows/naomi_ingest.py
do
    if [[ -e "$f" ]]; then
        run "git rm '$f' || rm -f '$f'"
    fi
done

# ------------------------------------------------------------------
# Batch 2: delete dashboard notebook + stale marimo caches
# ------------------------------------------------------------------

echo ""
echo "== Batch 2: delete dashboard + notebook caches =="
[[ -e notebooks/dashboard.py ]] && run "git rm notebooks/dashboard.py || rm -f notebooks/dashboard.py"
[[ -d notebooks/__marimo__ ]] && run "rm -rf notebooks/__marimo__"
[[ -d notebooks/__pycache__ ]] && run "rm -rf notebooks/__pycache__"

# ------------------------------------------------------------------
# Batch 3: delete out-of-scope data outputs
# ------------------------------------------------------------------

echo ""
echo "== Batch 3: delete out-of-scope data outputs =="
for d in data/gold data/models; do
    [[ -e "$d" ]] && run "rm -rf '$d'"
done
[[ -e data/interim/naomi_labels.parquet ]] && run "rm -f data/interim/naomi_labels.parquet"

# ------------------------------------------------------------------
# Batch 4: move causal/ to repo root (per CONTEXT.md monorepo layout)
# ------------------------------------------------------------------

echo ""
echo "== Batch 4: move causal/ to repo root =="
if [[ -d src/hivdrwatch/causal ]]; then
    run "git mv src/hivdrwatch/causal causal || mv src/hivdrwatch/causal causal"
    # Also move flows/causal.py alongside — it belongs with the causal package.
    if [[ -e src/hivdrwatch/flows/causal.py ]]; then
        run "mkdir -p causal/flows"
        run "git mv src/hivdrwatch/flows/causal.py causal/flows/causal.py \
             || mv src/hivdrwatch/flows/causal.py causal/flows/causal.py"
        run "touch causal/flows/__init__.py"
    fi
fi

# Move tests too
if [[ -d tests/causal ]]; then
    # Leave in place — they still work with updated imports.
    :
fi

# ------------------------------------------------------------------
# Batch 5: rename hivdb/ → resistance/, curate/ → align/, flows/ → pipeline/
# ------------------------------------------------------------------

echo ""
echo "== Batch 5: rename dirs to MODULES.md contract =="

rename_dir() {
    local src="$1"; local dst="$2"
    if [[ -d "$src" ]] && [[ ! -d "$dst" ]]; then
        run "git mv '$src' '$dst' || mv '$src' '$dst'"
    fi
}

rename_dir src/hivdrwatch/hivdb    src/hivdrwatch/resistance
rename_dir src/hivdrwatch/curate   src/hivdrwatch/align
rename_dir src/hivdrwatch/flows    src/hivdrwatch/pipeline

# Test dirs
rename_dir tests/hivdb   tests/resistance
rename_dir tests/curate  tests/align
rename_dir tests/flows   tests/pipeline

# ------------------------------------------------------------------
# Batch 6: create empty workflows/ for future Nextflow paper analyses
# ------------------------------------------------------------------

echo ""
echo "== Batch 6: scaffold workflows/ dir =="
if [[ ! -d workflows ]]; then
    run "mkdir -p workflows"
    if [[ $APPLY -eq 1 ]]; then
        cat > workflows/README.md <<'MD'
# workflows/

Nextflow (nf-core convention) pipelines for **paper analyses** — the
reproducible, one-command scripts that regenerate every figure from raw
data. Distinct from `src/hivdrwatch/pipeline/`, which is the live
surveillance service (Prefect 3).

See CLAUDE.md → "Choosing an orchestrator" for the rule.

Populated per paper as it's written:
  * `paper1_tools.nf` — full paper 1 pipeline (ingest → cluster → figures)
  * `paper2_causal.nf` — Paper 2 causal-identification pipeline
  * `paper3_foundation.nf` — Paper 3 benchmark pipeline
MD
    fi
fi

# ------------------------------------------------------------------
# Batch 7: rewrite imports across the repo
# ------------------------------------------------------------------

echo ""
echo "== Batch 7: rewrite imports for renames =="

# Substitutions to apply. Order matters — do longer patterns first.
declare -a SUBS=(
    's|from hivdrwatch\.hivdb|from hivdrwatch.resistance|g'
    's|import hivdrwatch\.hivdb|import hivdrwatch.resistance|g'
    's|from hivdrwatch\.curate|from hivdrwatch.align|g'
    's|import hivdrwatch\.curate|import hivdrwatch.align|g'
    's|from hivdrwatch\.flows|from hivdrwatch.pipeline|g'
    's|import hivdrwatch\.flows|import hivdrwatch.pipeline|g'
    's|from hivdrwatch\.causal|from causal|g'
    's|import hivdrwatch\.causal|import causal|g'
    # Relative imports inside the package: ..curate -> ..align, etc.
    's|from \.\.hivdb|from ..resistance|g'
    's|from \.\.curate|from ..align|g'
    's|from \.\.flows|from ..pipeline|g'
    # `hivdrwatch.flows.causal` inside modules that reference the CLI path
    's|hivdrwatch\.flows\.causal|causal.flows.causal|g'
)

# Find target files
FILES=$(find src tests scripts causal -type f \( -name "*.py" -o -name "*.md" \) 2>/dev/null | grep -v __pycache__ || true)

for sub in "${SUBS[@]}"; do
    for f in $FILES; do
        if [[ -f "$f" ]] && grep -qE "hivdrwatch\.(hivdb|curate|flows|causal)|\.\.(hivdb|curate|flows)" "$f" 2>/dev/null; then
            run "sed -i.bak '$sub' '$f' && rm '$f.bak'"
        fi
    done
done

# Update Makefile references to renamed CLI modules
if [[ -e Makefile ]]; then
    run "sed -i.bak 's|hivdrwatch\.flows\.|hivdrwatch.pipeline.|g' Makefile && rm Makefile.bak"
    # Remove Makefile targets that referenced deleted modules
    if [[ $APPLY -eq 1 ]]; then
        python3 - <<'PY'
import re, pathlib
mk = pathlib.Path("Makefile")
text = mk.read_text()
# Drop naomi, dr-features, train-dr, train-outbreak, dashboard targets.
for target in ("naomi", "dr-features", "train-dr", "train-outbreak", "dashboard", "ingest-full"):
    # Match `.PHONY: <target>\n<target>:...\n\t...` blocks including the recipe lines
    text = re.sub(
        rf"(?ms)^\.PHONY: {re.escape(target)}\n{re.escape(target)}:.*?(?=\n\n|\n\.PHONY|\Z)",
        "",
        text,
    )
mk.write_text(text)
print("stripped removed-module targets from Makefile")
PY
    fi
fi

# ------------------------------------------------------------------
# Batch 8: pyproject.toml + entry points
# ------------------------------------------------------------------

echo ""
echo "== Batch 8: update pyproject.toml =="
if [[ -e pyproject.toml ]] && [[ $APPLY -eq 1 ]]; then
    python3 - <<'PY'
import pathlib, re
p = pathlib.Path("pyproject.toml")
text = p.read_text()

# Update entry point that referenced hivdrwatch.hivdb.__main__
text = text.replace("hivdrwatch.hivdb.__main__", "hivdrwatch.resistance.__main__")

# Add causal/ as a package alongside src/hivdrwatch/
if 'packages = ["src/hivdrwatch"]' in text:
    text = text.replace(
        'packages = ["src/hivdrwatch"]',
        'packages = ["src/hivdrwatch", "causal"]',
    )

# Add "." to pytest pythonpath so top-level causal/ imports resolve
text = re.sub(
    r'pythonpath = \["src"\]',
    'pythonpath = ["src", "."]',
    text,
)

p.write_text(text)
print("pyproject.toml updated (entry point + packages + pytest pythonpath)")
PY
fi

# ------------------------------------------------------------------
# Done
# ------------------------------------------------------------------

echo ""
if [[ $APPLY -eq 1 ]]; then
    echo "== DONE =="
    echo ""
    echo "Next steps:"
    echo "  1. git status              # review the changes"
    echo "  2. pytest tests/ -q        # make sure imports still resolve"
    echo "  3. git add -A && git commit -m 'refactor: phase 0 — align repo to MODULES.md contract'"
    echo ""
    echo "If pytest fails on import errors, run:"
    echo "  grep -rn 'hivdrwatch\.(hivdb|curate|flows)' src tests"
    echo "and fix the leftover references."
else
    echo "== DRY-RUN COMPLETE — re-run with --apply to execute =="
fi
