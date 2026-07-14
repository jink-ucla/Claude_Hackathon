# DATA DICTIONARY

Column reference for every output CSV. Machine-readable shapes: `data/MANIFEST.json`.
Source legend: **DE** = `DE_stats.suppl_table.csv`; **KD** = `guide_kd_efficiency`; **REP** =
donor/guide correlation tables; **CLU** = `clustering_results_and_annotations`; **AUTO** =
`cluster_autoimmune_enrichment`; **POL/AGE** = polarization/aging coefficient tables;
**K562** = `K562_comparison`; **VAL** = `IL10_IL21_arrayed_validation`; **derived** = computed here.

---

## `data/integrated/regulator_master.csv` — ★ master (33,983 rows × 46 cols)

One row per **perturbed gene × culture condition**.

| Column | Type | Source | Description |
|---|---|---|---|
| `gene`, `gene_id` | str | DE | perturbed-gene symbol / Ensembl ID |
| `condition` | str | DE | Rest / Stim8hr / Stim48hr |
| `n_cells_target` | int | DE | cells carrying a targeting guide |
| `target_baseMean` | float | DE | mean baseline expression of the target |
| `ontarget_significant` | bool | DE | on-target knockdown significant (10% FDR) |
| `ontarget_effect_size` | float | DE | effect size on the intended target |
| `ontarget_effect_category` | str | DE | on-target KD / no on-target KD / putative off-target |
| `offtarget_flag` | bool | DE | potential off-target concern |
| `n_up_genes` / `n_down_genes` / `n_total_de_genes` | int | DE | significant DE genes (10% FDR) |
| `n_downstream` | int | DE | **trans-effects** (DE genes excluding on-target) — breadth metric |
| `n_total_genes_category` | str | DE | binned trans-effect size (`no effect` … `>10 DE genes`) |
| `is_strong_regulator` | bool | derived | `n_total_genes_category == ">10 DE genes"` |
| `n_guides_tested` / `n_guides_signif_kd` | int | KD | guides for this gene×condition / with sig. KD |
| `frac_guides_signif_kd` | float | derived | fraction of guides with significant knockdown |
| `best_kd_tstat` | float | KD | strongest (most negative) knockdown t-statistic |
| `is_efficient_kd` | bool | derived | `frac_guides_signif_kd ≥ 0.5` |
| `donor_correlation_mean` / `_min` | float | REP | cross-donor DE z-score correlation (NaN if untested) |
| `is_reproducible_donor` | bool | derived | `donor_correlation_mean ≥ 0.4` |
| `guide_correlation` | float | REP | cross-guide DE z-score correlation (NaN if single-guide) |
| `is_reproducible_guide` | bool | derived | `guide_correlation ≥ 0.4` |
| `cluster` | float | CLU | HDBSCAN GRN module id (NaN if not a module member) |
| `cluster_annotation` | str | CLU | manual module annotation (e.g. "Th17 differentiation") |
| `condition_specificity` | str | CLU | module's condition specificity |
| `in_grn_module` | bool | derived | gene×condition assigned to a GRN module |
| `autoimmune_or` / `autoimmune_fdr` | float | AUTO | module's autoimmune-disease odds ratio / BH-FDR |
| `is_autoimmune_enriched` | bool | derived | `autoimmune_or > 1 & autoimmune_fdr < 0.1` |
| `activation_coef_rank` / `polarization_coef_rank` | float | POL | activation / Th2-vs-Th1 regulator rank (0–1) |
| `is_known_polarization_reg` | bool | POL | annotated known activation/polarization regulator |
| `is_polarization_regulator` | bool | derived | `coef_rank ≥ 0.9` (either signature) or known |
| `aging_coef_rank` | float | AGE | CD4T aging-signature regulator rank (0–1) |
| `is_known_aging_reg` | bool | AGE | annotated known aging regulator |
| `is_aging_regulator` | bool | derived | `aging_coef_rank ≥ 0.9` or known |
| `k562_logfc_r` | float | K562 | K562-vs-CD4 logFC Pearson r (NaN if untested) |
| `n_degs_mash` | float | K562 | MASH DEG count in this CD4 condition |
| `is_cd4_specific` | bool | derived | `k562_logfc_r < 0.2`; NaN where untested (K562 comparison absent) |
| `il10_log2fc_vs_ntc` / `il21_log2fc_vs_ntc` | float | VAL | arrayed IL10/IL21 log2FC vs NTC (NaN if not validated) |
| `is_il_validated` | bool | derived | perturbation present in arrayed validation |
| `regulator_score` | int | derived | convergent-evidence score (points table below) |

> **`regulator_score` (0–9) — additive, transparent, reweightable:**
> `is_strong_regulator ×2` + `ontarget_significant` + `is_efficient_kd` + `is_reproducible_donor`
> + `is_reproducible_guide` + `is_polarization_regulator` + `is_aging_regulator`
> + `is_autoimmune_enriched`. A **ranking heuristic, not a statistical model** — each term is one
> independent evidence layer, so the score rewards convergence; recompute from the columns above
> to reweight. Guaranteed: `is_strong_regulator ⇒ score ≥ 2`; `ontarget_significant ⇒ score ≥ 1`.

## `data/integrated/condition_summary.csv` (3 rows)
Per-condition rollup: `n_perturbations`, `n_ontarget_significant`, `n_strong_regulator`,
`median_n_downstream`, `n_efficient_kd`, `n_reproducible_donor`, `n_polarization_reg`,
`n_aging_reg`, `n_autoimmune_enriched`, `n_score_ge_5`, `n_score_ge_6`.

## `data/integrated/cluster_summary.csv` (112 rows)
Per GRN module: `cluster`, `manual_annotation`, `best_described_by`, `condition_specificity`,
`cluster_size`, membership rollup (`n_member_genes`, `mean_regulator_score`, `top_members`),
`best_autoimmune_or`, `min_autoimmune_fdr`, `is_autoimmune_enriched`.

## `data/clean/*.csv`
Tidy one-table-per-source intermediates (see `data/MANIFEST.json` for each schema):
`de_stats`, `kd_efficiency_by_gene`, `donor_reproducibility`, `guide_reproducibility`,
`cluster_annotations`, `cluster_membership`, `cluster_autoimmune`, `polarization_regulators`,
`aging_regulators`, `k562_comparison`, `il10il21_validation`.
