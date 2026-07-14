#!/usr/bin/env python
"""
build_dataset.py — CD4+ T cell Perturb-seq master-regulator integration.

Integrates the shipped functional-genomics layers of the Marson/Pritchard 2025
genome-scale CRISPRi Perturb-seq screen in primary human CD4+ T cells into ONE
scored master table, keyed by (perturbed gene x culture condition) — because the
paper's central thesis is that active regulators change across stimulation
conditions, so the table is stratified by Rest / Stim8hr / Stim48hr, never pooled.

Evidence layers (each an independent axis; the score rewards convergence):
  - DE / trans-effect breadth   (DE_stats.suppl_table.csv)
  - CRISPRi knockdown efficiency (guide_kd_efficiency.suppl_table.csv)
  - cross-donor reproducibility  (DE_donor_robustness_correlation_summary.csv)
  - cross-guide reproducibility  (DE_by_guide.correlation_results.csv)
  - GRN module membership        (clustering_results_and_annotations.csv)
  - autoimmune-disease enrichment(cluster_autoimmune_enrichment_results.suppl_table.csv)
  - polarization/activation role (polarization_prediction_..._coefficients.csv)
  - aging role                   (aging_prediction_..._coefficients.csv)
  - K562 cross-cell-type spec.   (K562_comparison.suppl_table.csv)
  - arrayed IL10/IL21 validation (IL10_IL21_arrayed_validation.csv)

Run: <venv>/bin/python scripts/build_dataset.py
"""
from __future__ import annotations
import os
import sys
import ast
import json
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------- paths
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (portable)
RAW = os.path.join(ROOT, "suppl_tables")
CLEAN = os.path.join(ROOT, "data", "clean")
INTEG = os.path.join(ROOT, "data", "integrated")
for d in (CLEAN, INTEG):
    os.makedirs(d, exist_ok=True)

CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]

# ---------------------------------------------------------------- thresholds (reweightable)
STRONG_CATEGORY = ">10 DE genes"   # DE_stats n_total_genes_category tier = strong regulator
KD_FRAC_MIN = 0.5                  # >= half of a gene's guides show significant knockdown
DONOR_CORR_MIN = 0.4               # cross-donor DE z-score correlation (mean over donor pairs)
GUIDE_CORR_MIN = 0.4               # cross-guide DE z-score correlation
COEF_RANK_MIN = 0.9               # regulator coef_rank (0-1) tier for polarization/aging
AUTOIMMUNE_OR_MIN = 1.0           # cluster autoimmune odds ratio
AUTOIMMUNE_FDR_MAX = 0.10         # cluster autoimmune BH-FDR
CD4_SPECIFIC_R_MAX = 0.2          # K562-vs-CD4 logFC Pearson r below this => CD4-specific

# ---------------------------------------------------------------- helpers
manifest: dict = {}


def rp(name: str) -> str:
    return os.path.join(RAW, name)


def _save(df: pd.DataFrame, relpath: str, note: str) -> pd.DataFrame:
    path = os.path.join(ROOT, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False)
    manifest[relpath] = {
        "rows": int(df.shape[0]),
        "cols": int(df.shape[1]),
        "columns": list(map(str, df.columns)),
        "note": note,
    }
    print(f"  [save] {relpath:52s} {df.shape[0]:>7d} x {df.shape[1]:<3d}  {note}")
    return df


def join(vals) -> str:
    xs = sorted({str(v) for v in vals if pd.notna(v) and str(v) != ""})
    return ";".join(xs)


def cond_gene_set(cond: str) -> str:
    return f"downstream_{cond}"


# ================================================================ [1] clean tables
def clean_de() -> pd.DataFrame:
    df = pd.read_csv(rp("DE_stats.suppl_table.csv"))
    df = df.rename(columns={"index": "obs_key",
                            "target_contrast_gene_name": "gene",
                            "target_contrast": "gene_id",
                            "culture_condition": "condition"})
    df["is_strong_regulator"] = df["n_total_genes_category"].eq(STRONG_CATEGORY)
    keep = ["obs_key", "gene", "gene_id", "condition", "n_cells_target", "target_baseMean",
            "n_up_genes", "n_down_genes", "n_total_de_genes", "n_downstream",
            "n_total_genes_category", "ontarget_effect_size", "ontarget_significant",
            "ontarget_effect_category", "offtarget_flag", "is_strong_regulator"]
    df = df[keep]
    return _save(df, "data/clean/de_stats.csv",
                 "one row per perturbation x condition; trans-effect breadth + on-target QC")


def clean_kd() -> pd.DataFrame:
    df = pd.read_csv(rp("guide_kd_efficiency.suppl_table.csv"))
    df = df.rename(columns={df.columns[0]: "guide_id"})
    agg = (df.groupby(["perturbed_gene_id", "culture_condition"], as_index=False)
             .agg(n_guides_tested=("signif_knockdown", "size"),
                  n_guides_signif_kd=("signif_knockdown", "sum"),
                  best_kd_tstat=("t_statistic", "min")))
    agg["frac_guides_signif_kd"] = agg["n_guides_signif_kd"] / agg["n_guides_tested"]
    agg = agg.rename(columns={"perturbed_gene_id": "gene_id", "culture_condition": "condition"})
    return _save(agg, "data/clean/kd_efficiency_by_gene.csv",
                 "CRISPRi knockdown efficiency aggregated to gene x condition")


def clean_donor() -> pd.DataFrame:
    df = pd.read_csv(rp("DE_donor_robustness_correlation_summary.csv"))
    df = df.rename(columns={"target_contrast": "gene_id", "target_name": "gene"})
    df = df[["gene", "gene_id", "condition", "donor_correlation_mean", "donor_correlation_min"]]
    df = df.drop_duplicates(["gene_id", "condition"])
    return _save(df, "data/clean/donor_reproducibility.csv",
                 "cross-donor DE z-score correlation per gene x condition")


def clean_guidecorr() -> pd.DataFrame:
    df = pd.read_csv(rp("DE_by_guide.correlation_results.csv"))
    df = df.rename(columns={"target": "gene", "culture_condition": "condition",
                            "correlation": "guide_correlation"})
    df = df[["gene", "condition", "guide_correlation", "correlation_ceiling", "n_signif_union"]]
    df = df.drop_duplicates(["gene", "condition"])
    return _save(df, "data/clean/guide_reproducibility.csv",
                 "cross-guide DE z-score correlation per gene x condition")


def clean_clusters() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(rp("clustering_results_and_annotations.csv"), encoding="utf-8-sig")
    df = df.rename(columns={df.columns[0]: "cluster"})
    ann_cols = ["cluster", "manual_annotation", "best_described_by", "condition_specificity",
                "cluster_size", "cluster_gene_size", "intracluster_corr",
                "rest_count", "stim8hr_count", "stim48hr_count"]
    ann = df[[c for c in ann_cols if c in df.columns]].copy()
    _save(ann, "data/clean/cluster_annotations.csv",
          "one row per HDBSCAN GRN cluster with annotation + condition specificity")

    # explode cluster_member_with_condition -> (gene, condition, cluster)
    recs = []
    for _, r in df.iterrows():
        try:
            members = ast.literal_eval(r["cluster_member_with_condition"])
        except Exception:
            continue
        for m in members:
            g, _, c = str(m).rpartition("_")
            if c in CONDITIONS and g:
                recs.append((g, c, int(r["cluster"])))
    mem = pd.DataFrame(recs, columns=["gene", "condition", "cluster"]).drop_duplicates()
    mem = mem.merge(ann[["cluster", "manual_annotation", "best_described_by",
                         "condition_specificity"]], on="cluster", how="left")
    # keep one cluster per (gene, condition) — HDBSCAN assigns a single module
    mem = mem.drop_duplicates(["gene", "condition"])
    _save(mem, "data/clean/cluster_membership.csv",
          "GRN module membership per gene x condition")
    return ann, mem


def clean_autoimmune() -> pd.DataFrame:
    df = pd.read_csv(rp("cluster_autoimmune_enrichment_results.suppl_table.csv"))
    df = df[(df["disease"] == "autoimmune disease") & (~df["negative_control_disease"].astype(bool))]
    df = df[df["gene_set"].isin([cond_gene_set(c) for c in CONDITIONS])].copy()
    df["condition"] = df["gene_set"].str.replace("downstream_", "", regex=False)
    df = df.rename(columns={"odds_ratio": "autoimmune_or", "p_adj_fdr": "autoimmune_fdr"})
    df = df[["cluster", "condition", "autoimmune_or", "autoimmune_fdr", "cluster_size",
             "in_cluster_in_disease", "intersecting_genes"]]
    df = df.drop_duplicates(["cluster", "condition"])
    return _save(df, "data/clean/cluster_autoimmune.csv",
                 "autoimmune-disease enrichment of each GRN cluster's downstream genes, per condition")


def clean_polarization() -> pd.DataFrame:
    df = pd.read_csv(rp("polarization_prediction_condition_comparison_regulator_coefficients.csv"))
    df = df.rename(columns={df.columns[0]: "row_key", "celltype": "condition"})
    df = df[df["condition"].isin(CONDITIONS)].copy()
    # signature: 'activation' and 'ota' (Th2-vs-Th1 polarization, Ota 2021)
    df["signature"] = df["signature"].replace({"ota": "polarization"})
    piv = df.pivot_table(index=["regulator", "condition"], columns="signature",
                         values="coef_rank", aggfunc="max")
    piv.columns = [f"{c}_coef_rank" for c in piv.columns]
    known = (df.groupby(["regulator", "condition"])["known_regulators"].any()
               .rename("is_known_polarization_reg"))
    out = piv.join(known).reset_index().rename(columns={"regulator": "gene"})
    return _save(out, "data/clean/polarization_regulators.csv",
                 "activation & Th2/Th1-polarization regulator coef_rank per gene x condition")


def clean_aging() -> pd.DataFrame:
    df = pd.read_csv(rp("aging_prediction_condition_comparison_regulator_coefficients.csv"))
    df = df.rename(columns={df.columns[0]: "row_key", "celltype": "condition"})
    df = df[df["condition"].isin(CONDITIONS)].copy()
    out = (df.groupby(["regulator", "condition"], as_index=False)
             .agg(aging_coef_rank=("coef_rank", "max"),
                  is_known_aging_reg=("known_regulators", "any"))
             .rename(columns={"regulator": "gene"}))
    return _save(out, "data/clean/aging_regulators.csv",
                 "CD4T aging-signature regulator coef_rank per gene x condition")


def clean_k562() -> pd.DataFrame:
    df = pd.read_csv(rp("K562_comparison.suppl_table.csv"))
    df = df.rename(columns={"target_contrast_gene_name": "gene"})
    ndeg = df.apply(lambda r: r.get(f"n_degs_MASH_{r['condition']}", np.nan), axis=1)
    df = df.assign(n_degs_mash=ndeg)
    df = df[["gene", "condition", "logfc_pearson_r", "donor_correlation_mean", "n_degs_mash"]]
    df = df.rename(columns={"logfc_pearson_r": "k562_logfc_r",
                            "donor_correlation_mean": "k562_donor_corr_mean"})
    df = df.drop_duplicates(["gene", "condition"])
    return _save(df, "data/clean/k562_comparison.csv",
                 "K562-vs-CD4 logFC correlation (cell-type specificity) per gene x condition")


def clean_validation() -> pd.DataFrame:
    df = pd.read_csv(rp("IL10_IL21_arrayed_validation.csv"), encoding="utf-8-sig")
    df = df.rename(columns={df.columns[0]: "sample"})
    ntc = df[df["Perturbation"] == "NTC"][["IL10_perc", "IL21_perc"]].mean()
    agg = (df.groupby("Perturbation", as_index=False)
             .agg(il10_perc=("IL10_perc", "mean"), il21_perc=("IL21_perc", "mean"),
                  n_donors_val=("Donor", "nunique")))
    agg["il10_log2fc_vs_ntc"] = np.log2((agg["il10_perc"] + 1e-3) / (ntc["IL10_perc"] + 1e-3))
    agg["il21_log2fc_vs_ntc"] = np.log2((agg["il21_perc"] + 1e-3) / (ntc["IL21_perc"] + 1e-3))
    agg = agg.rename(columns={"Perturbation": "gene"})
    return _save(agg, "data/clean/il10il21_validation.csv",
                 "arrayed flow-cytometry IL10/IL21 validation per perturbation (all conditions)")


# ================================================================ [2] integrate
def build_master(de, kd, donor, guidecorr, mem, auto, pol, aging, k562, val) -> pd.DataFrame:
    m = de.copy()

    # --- KD efficiency (gene_id x condition)
    m = m.merge(kd, on=["gene_id", "condition"], how="left")
    m["frac_guides_signif_kd"] = m["frac_guides_signif_kd"].fillna(0.0)
    m["is_efficient_kd"] = m["frac_guides_signif_kd"] >= KD_FRAC_MIN

    # --- reproducibility
    m = m.merge(donor[["gene_id", "condition", "donor_correlation_mean", "donor_correlation_min"]],
                on=["gene_id", "condition"], how="left")
    m["is_reproducible_donor"] = m["donor_correlation_mean"] >= DONOR_CORR_MIN
    m = m.merge(guidecorr[["gene", "condition", "guide_correlation"]],
                on=["gene", "condition"], how="left")
    m["is_reproducible_guide"] = m["guide_correlation"] >= GUIDE_CORR_MIN

    # --- GRN module membership + autoimmune enrichment via cluster
    m = m.merge(mem[["gene", "condition", "cluster", "manual_annotation",
                     "condition_specificity"]], on=["gene", "condition"], how="left")
    m = m.rename(columns={"manual_annotation": "cluster_annotation"})
    m["in_grn_module"] = m["cluster"].notna()
    m = m.merge(auto[["cluster", "condition", "autoimmune_or", "autoimmune_fdr"]],
                on=["cluster", "condition"], how="left")
    m["is_autoimmune_enriched"] = ((m["autoimmune_or"] > AUTOIMMUNE_OR_MIN) &
                                   (m["autoimmune_fdr"] < AUTOIMMUNE_FDR_MAX)).fillna(False)

    # --- polarization / activation regulator
    m = m.merge(pol, on=["gene", "condition"], how="left")
    for c in ["activation_coef_rank", "polarization_coef_rank"]:
        if c not in m.columns:
            m[c] = np.nan
    m["is_known_polarization_reg"] = m["is_known_polarization_reg"].fillna(False)
    m["is_polarization_regulator"] = (
        (m[["activation_coef_rank", "polarization_coef_rank"]].max(axis=1) >= COEF_RANK_MIN)
        | m["is_known_polarization_reg"]).fillna(False)

    # --- aging regulator
    m = m.merge(aging, on=["gene", "condition"], how="left")
    m["is_known_aging_reg"] = m["is_known_aging_reg"].fillna(False)
    m["is_aging_regulator"] = ((m["aging_coef_rank"] >= COEF_RANK_MIN)
                               | m["is_known_aging_reg"]).fillna(False)

    # --- K562 cross-cell-type specificity
    m = m.merge(k562, on=["gene", "condition"], how="left")
    m["is_cd4_specific"] = m["k562_logfc_r"].lt(CD4_SPECIFIC_R_MAX).where(m["k562_logfc_r"].notna())  # NaN where untested

    # --- arrayed IL10/IL21 validation (per gene, all conditions)
    m = m.merge(val[["gene", "il10_log2fc_vs_ntc", "il21_log2fc_vs_ntc", "n_donors_val"]],
                on="gene", how="left")
    m["is_il_validated"] = m["n_donors_val"].notna()

    # --- convergent evidence score (0-9)
    b = lambda s: s.fillna(False).astype(int)
    m["regulator_score"] = (
        b(m["is_strong_regulator"]) * 2
        + b(m["ontarget_significant"])
        + b(m["is_efficient_kd"])
        + b(m["is_reproducible_donor"])
        + b(m["is_reproducible_guide"])
        + b(m["is_polarization_regulator"])
        + b(m["is_aging_regulator"])
        + b(m["is_autoimmune_enriched"])
    )

    cols = ["gene", "gene_id", "condition",
            "n_cells_target", "target_baseMean", "ontarget_significant",
            "ontarget_effect_size", "ontarget_effect_category", "offtarget_flag",
            "n_up_genes", "n_down_genes", "n_total_de_genes", "n_downstream",
            "n_total_genes_category", "is_strong_regulator",
            "n_guides_tested", "n_guides_signif_kd", "frac_guides_signif_kd",
            "best_kd_tstat", "is_efficient_kd",
            "donor_correlation_mean", "donor_correlation_min", "is_reproducible_donor",
            "guide_correlation", "is_reproducible_guide",
            "cluster", "cluster_annotation", "condition_specificity", "in_grn_module",
            "autoimmune_or", "autoimmune_fdr", "is_autoimmune_enriched",
            "activation_coef_rank", "polarization_coef_rank",
            "is_known_polarization_reg", "is_polarization_regulator",
            "aging_coef_rank", "is_known_aging_reg", "is_aging_regulator",
            "k562_logfc_r", "n_degs_mash", "is_cd4_specific",
            "il10_log2fc_vs_ntc", "il21_log2fc_vs_ntc", "is_il_validated",
            "regulator_score"]
    m = m[cols].sort_values(["regulator_score", "n_downstream"], ascending=[False, False])
    return _save(m, "data/integrated/regulator_master.csv",
                 "MASTER: one row per perturbed gene x condition, all evidence layers + regulator_score")


def build_condition_summary(m) -> pd.DataFrame:
    rows = []
    for c in CONDITIONS:
        s = m[m["condition"] == c]
        rows.append(dict(
            condition=c,
            n_perturbations=len(s),
            n_ontarget_significant=int(s["ontarget_significant"].sum()),
            n_strong_regulator=int(s["is_strong_regulator"].sum()),
            median_n_downstream=float(s["n_downstream"].median()),
            n_efficient_kd=int(s["is_efficient_kd"].sum()),
            n_reproducible_donor=int(s["is_reproducible_donor"].sum()),
            n_polarization_reg=int(s["is_polarization_regulator"].sum()),
            n_aging_reg=int(s["is_aging_regulator"].sum()),
            n_autoimmune_enriched=int(s["is_autoimmune_enriched"].sum()),
            n_score_ge_5=int((s["regulator_score"] >= 5).sum()),
            n_score_ge_6=int((s["regulator_score"] >= 6).sum()),
        ))
    return _save(pd.DataFrame(rows), "data/integrated/condition_summary.csv",
                 "per-condition rollup of the master table")


def build_cluster_summary(m, ann, auto) -> pd.DataFrame:
    memc = m[m["in_grn_module"]].groupby("cluster").agg(
        n_member_gene_conditions=("gene", "size"),
        n_member_genes=("gene", "nunique"),
        mean_regulator_score=("regulator_score", "mean"),
        top_members=("gene", lambda s: join(list(pd.Series(s).drop_duplicates())[:8])),
    )
    au = (auto.sort_values("autoimmune_fdr")
              .groupby("cluster")
              .agg(best_autoimmune_or=("autoimmune_or", "max"),
                   min_autoimmune_fdr=("autoimmune_fdr", "min")))
    out = (ann.set_index("cluster")
              .join(memc).join(au).reset_index())
    out["is_autoimmune_enriched"] = ((out["best_autoimmune_or"] > AUTOIMMUNE_OR_MIN) &
                                     (out["min_autoimmune_fdr"] < AUTOIMMUNE_FDR_MAX)).fillna(False)
    out = out.sort_values("min_autoimmune_fdr")
    return _save(out, "data/integrated/cluster_summary.csv",
                 "per GRN cluster: annotation, membership, autoimmune enrichment")


# ================================================================ main
def main():
    print("[1] clean tables")
    de = clean_de()
    kd = clean_kd()
    donor = clean_donor()
    guidecorr = clean_guidecorr()
    ann, mem = clean_clusters()
    auto = clean_autoimmune()
    pol = clean_polarization()
    aging = clean_aging()
    k562 = clean_k562()
    val = clean_validation()

    print("[2] integrated tables")
    master = build_master(de, kd, donor, guidecorr, mem, auto, pol, aging, k562, val)
    build_condition_summary(master)
    build_cluster_summary(master, ann, auto)

    with open(os.path.join(ROOT, "data", "MANIFEST.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  [save] data/MANIFEST.json  ({len(manifest)} tables)")

    # cheap invariants
    assert master.shape[0] == de.shape[0] == 33983, f"master rows {master.shape[0]}"
    assert not master.duplicated(["gene_id", "condition"]).any(), "duplicate gene_id x condition"
    assert master["regulator_score"].between(0, 9).all(), "score out of [0,9]"
    print(f"\n[ok] master = {master.shape[0]} rows; "
          f"score>=6: {(master.regulator_score>=6).sum()}; "
          f"top gene(Stim8hr): {master[master.condition=='Stim8hr'].iloc[0]['gene']}")


if __name__ == "__main__":
    sys.exit(main())
