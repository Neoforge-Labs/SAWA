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
