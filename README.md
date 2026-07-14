# CD4⁺ T cell Regulator Atlas — genome-scale Perturb-seq integration

**One scored, honestly-ranked master table of which genes rewire the human CD4⁺ T cell
transcriptome when knocked down — stratified by stimulation state and linked to helper-cell
polarization, aging, and autoimmune-disease genetics.**

Built from the processed tables of the Marson/Pritchard 2025 genome-scale CRISPRi Perturb-seq
screen in primary human CD4⁺ T cells (~22M cells, 4 donors, Rest/Stim8hr/Stim48hr;
bioRxiv [10.64898/2025.12.23.696273](https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1),
data `emdann/GWT_perturbseq_analysis_2025`). No raw single cells required — the integration,
prioritization, and honest scoring are the value-add.

## What this is — eight evidence layers, one table

| Layer | What's in the data | Source table |
|---|---|---|
| DE / trans-effect breadth | up/down/total DE genes, `n_downstream` | `DE_stats.suppl_table.csv` |
| On-target validity | `ontarget_significant`, KD efficiency | `DE_stats`, `guide_kd_efficiency` |
| Reproducibility | cross-guide & cross-donor DE correlation | `DE_by_guide.*`, `DE_donor_robustness_*` |
| GRN module | HDBSCAN cluster + annotation | `clustering_results_and_annotations.csv` |
| Autoimmune enrichment | per-module odds ratio + FDR | `cluster_autoimmune_enrichment_results.*` |
| Polarization/activation | regulator `coef_rank` | `polarization_prediction_*_coefficients.csv` |
| Aging | CD4T aging regulator `coef_rank` | `aging_prediction_*_coefficients.csv` |
| K562 specificity | K562-vs-CD4 logFC correlation | `K562_comparison.suppl_table.csv` |

Provenance and exact thresholds: `docs/GROUNDING.md`. Every column: `docs/DATA_DICTIONARY.md`.

## Layout

```
Claude_Hackathon/
├── analysis.md              # findings narrative (start here)
├── README.md
├── SUBMISSION.md            # hackathon submission write-up / abstract
├── LICENSE                  # MIT (+ upstream data attribution)
├── data_sharing_readme.md   # upstream data-sharing notes (Marson/Pritchard)
├── suppl_tables/            # raw inputs (shipped supplementary CSVs)
├── metadata/                # per-file provenance (12 *.jsonld)
├── data/
│   ├── clean/               # 11 tidy per-source CSVs
│   ├── integrated/
│   │   ├── regulator_master.csv   # ★ 33,983 gene×condition × 46 cols (incl. regulator_score)
│   │   ├── condition_summary.csv  # 3-row per-condition rollup
│   │   ├── cluster_summary.csv    # 112 GRN modules + autoimmune enrichment
│   │   └── perturbation_baseline_summary.csv  # cz-benchmarks-idiom baseline scores (Task 2)
│   └── MANIFEST.json        # machine-readable schema of every output
├── figures/                 # 7 PNGs + regulator_module_network.graphml
├── scripts/                 # build_dataset · make_figures · validate_dataset · perturbation_baseline
└── docs/                    # GROUNDING.md · DATA_DICTIONARY.md · CD4_perturbseq_research_dossier.md
```

## The star table — `data/integrated/regulator_master.csv`

One row per **perturbed gene × culture condition**. Key columns: `gene`, `condition`,
`n_downstream` (trans-effect breadth), `ontarget_significant`, `frac_guides_signif_kd`,
`donor_correlation_mean`, `guide_correlation`, `cluster_annotation`, `autoimmune_or`,
`polarization_coef_rank`, `aging_coef_rank`, `k562_logfc_r`, and the additive **`regulator_score`
(0–9)**. Sorting by `regulator_score` surfaces top-scoring, data-nominated CD4⁺ regulator
candidates first — **BATF**, a canonical helper-cell master TF, and **BAHD1**, a data-nominated
chromatin factor (MLL2 complex), both at 9; STAT3, STAT6, GATA3, RASA2 at 8; the proximal TCR
machinery (CD3E, LAT, ZAP70, PLCG1) as the broadest trans-effect hubs — purely from data
integration. See `analysis.md` §2–3.

**Honesty notes (see analysis.md §2, §9):** lineage-definers FOXP3 (2) and RORC (1) score low —
these culture conditions capture activation/cytokine circuitry, not Treg/Th17 commitment; the
score reflects *evidence in this dataset*, not ground-truth importance.

## Reproduce

```bash
VP=./.venv312/bin/python                 # python3.12 venv: pandas numpy matplotlib networkx scipy anndata h5py
$VP scripts/build_dataset.py             # raw CSVs -> data/clean -> data/integrated -> MANIFEST.json
$VP scripts/make_figures.py              # -> figures/*.png + .graphml
$VP scripts/validate_dataset.py          # 36/36 integrity checks (raw vs built)
$VP scripts/perturbation_baseline.py     # cz-benchmarks-idiom prediction baseline (Task 2)
```

## Provenance & citation

Zhu, Dann, Yan, Reyes Retana, Goto, Guitche, Petersen, Ota, Pritchard, Marson.
*Genome-scale perturb-seq in primary human CD4⁺ T cells maps context-specific regulators of T
cell programs and human immune traits.* bioRxiv 2025. MIT License. Supplementary tables from
GitHub `emdann/GWT_perturbseq_analysis_2025`. Companion research dossier:
`docs/CD4_perturbseq_research_dossier.md`.
