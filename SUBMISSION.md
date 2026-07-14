# Submission — CD4⁺ T cell Regulator Atlas (genome-scale Perturb-seq integration)

**One-line:** an integrated, honestly-scored prioritization of gene-regulatory "master
regulators" in primary human CD4⁺ T cells, built from the Marson/Pritchard 2025 genome-scale
CRISPRi Perturb-seq screen (CZI Virtual Cells Platform dataset), stratified by stimulation state
and linked to helper-cell polarization, aging, and autoimmune-disease genetics.

## What it is
Eight functional-genomics evidence layers from the shipped screen are joined into one master
table — **33,983 perturbed-gene × condition rows** — with a transparent additive
**`regulator_score` (0–9)** that rewards convergent evidence (trans-effect breadth, on-target
QC, cross-guide/cross-donor reproducibility, GRN-module membership, autoimmune enrichment,
polarization/aging regulator role, K562 cross-cell-type specificity). A companion module
reimplements the cz-benchmarks *Perturbation Expression Prediction* metric as a model-free
baseline.

## Key results
- **The ranking recovers known biology unsupervised:** top regulators are BATF (9), STAT3 /
  STAT6 / GATA3 / RASA2 (8); the broadest trans-effect hubs are the proximal TCR machinery
  (CD3E, LAT, ZAP70, PLCG1, >5,000 downstream genes each).
- **Context-specificity, quantified:** the same perturbation's effect transfers across
  stimulation states only modestly (median Spearman ρ ≈ 0.20, ~20× the random baseline) and
  degrades with context distance — a direct confirmation of the paper's thesis.
- **Disease link:** autoimmune-enriched GRN modules are most numerous at Stim8hr (13) and
  Stim48hr (9), several of them novel ("unknown"-annotated) programs.
- **Rigor:** `validate_dataset.py` runs 36 raw-vs-built integrity checks (all pass); every
  score-feeding layer is re-derived from the raw tables.

## What's in this bundle
```
analysis.md                       ← findings narrative (read first)
README.md                         ← framing, layout, reproduce
SUBMISSION.md                     ← this file
data/integrated/regulator_master.csv        ← ★ the result (33,983 × 46 cols, incl. regulator_score)
data/integrated/{condition,cluster,perturbation_baseline}_summary.csv
data/clean/*.csv                  ← 11 tidy per-source tables
data/MANIFEST.json                ← machine-readable schema of every output
figures/*.png (7) + regulator_module_network.graphml
scripts/{build_dataset,make_figures,validate_dataset,perturbation_baseline}.py
docs/{GROUNDING,DATA_DICTIONARY,CD4_perturbseq_research_dossier}.md
suppl_tables/                     ← shipped input CSVs (inputs to the pipeline)
```
Not included (fetch on demand): `GWCD4i.DE_stats.h5ad` (16 GB, only needed for the Task-2
baseline) via `aws s3 cp --no-sign-request s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad .`
and the local Python venv.

## Reproduce
```bash
python3 -m venv .venv && .venv/bin/pip install pandas numpy matplotlib networkx scipy anndata h5py awscli
.venv/bin/python scripts/build_dataset.py       # -> data/clean, data/integrated, MANIFEST.json
.venv/bin/python scripts/make_figures.py        # -> figures/
.venv/bin/python scripts/validate_dataset.py    # 36/36 integrity checks
.venv/bin/python scripts/perturbation_baseline.py   # needs GWCD4i.DE_stats.h5ad
```

## Data & attribution
Zhu, Dann, Yan, Reyes Retana, Goto, Guitche, Petersen, Ota, Pritchard, Marson. *Genome-scale
perturb-seq in primary human CD4⁺ T cells maps context-specific regulators of T cell programs
and human immune traits.* bioRxiv 2025 (DOI 10.64898/2025.12.23.696273). MIT License.
Supplementary tables: GitHub `emdann/GWT_perturbseq_analysis_2025`. Caveats and honesty notes:
`analysis.md` §9 and `docs/GROUNDING.md` §5.
