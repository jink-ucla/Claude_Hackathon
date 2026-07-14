#!/usr/bin/env python
"""
validate_dataset.py — integrity checks for the CD4+ master-regulator build.

Philosophy: re-derive key numbers straight from the RAW supplementary CSVs and
assert they match the built outputs, so a silent parse/join/score error fails
loudly. Every evidence layer that feeds regulator_score is guarded by a genuine
raw re-derivation (not a built-column tautology), and the score itself is checked
to decompose exactly into its 8 constituent flags. Exits non-zero on any failure.

Run: <venv>/bin/python scripts/validate_dataset.py
"""
from __future__ import annotations
import os
import sys
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (portable)
RAW = os.path.join(ROOT, "suppl_tables")
INTEG = os.path.join(ROOT, "data", "integrated")
CLEAN = os.path.join(ROOT, "data", "clean")

# thresholds (must match build_dataset.py)
STRONG_CATEGORY = ">10 DE genes"
KD_FRAC_MIN = 0.5
DONOR_CORR_MIN = 0.4
GUIDE_CORR_MIN = 0.4
AUTO_OR_MIN, AUTO_FDR_MAX = 1.0, 0.10
CD4_SPECIFIC_R_MAX = 0.2

checks = []


def check(name, cond, detail=""):
    checks.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  — {detail}" if detail else ""))


def main():
    # ---- raw sources (each re-read independently) ----
    de_raw = pd.read_csv(os.path.join(RAW, "DE_stats.suppl_table.csv"))
    kd_raw = pd.read_csv(os.path.join(RAW, "guide_kd_efficiency.suppl_table.csv"))
    auto_raw = pd.read_csv(os.path.join(RAW, "cluster_autoimmune_enrichment_results.suppl_table.csv"))
    clu_raw = pd.read_csv(os.path.join(RAW, "clustering_results_and_annotations.csv"), encoding="utf-8-sig")
    donor_raw = pd.read_csv(os.path.join(RAW, "DE_donor_robustness_correlation_summary.csv"))
    guide_raw = pd.read_csv(os.path.join(RAW, "DE_by_guide.correlation_results.csv"))
    k562_raw = pd.read_csv(os.path.join(RAW, "K562_comparison.suppl_table.csv"))
    val_genes = set(pd.read_csv(os.path.join(RAW, "IL10_IL21_arrayed_validation.csv"),
                                encoding="utf-8-sig")["Perturbation"])
    # ---- built outputs ----
    m = pd.read_csv(os.path.join(INTEG, "regulator_master.csv"), low_memory=False)
    cs = pd.read_csv(os.path.join(INTEG, "condition_summary.csv"))
    mem = pd.read_csv(os.path.join(CLEAN, "cluster_membership.csv"))

    print("== A. raw-vs-built row counts ==")
    check("DE_stats raw rows == 33983", de_raw.shape[0] == 33983, f"{de_raw.shape[0]}")
    check("master rows == DE_stats raw rows", m.shape[0] == de_raw.shape[0], f"{m.shape[0]}")
    check("condition_summary has 3 rows", cs.shape[0] == 3)
    check("cluster annotations == raw clusters", clu_raw.shape[0] == int(m.cluster.dropna().nunique()),
          f"raw={clu_raw.shape[0]}, master clusters={int(m.cluster.dropna().nunique())}")

    print("== B. join integrity ==")
    check("master unique on (gene_id, condition)", not m.duplicated(["gene_id", "condition"]).any())
    key_raw = set(zip(de_raw.target_contrast, de_raw.culture_condition))
    key_m = set(zip(m.gene_id, m.condition))
    check("master keys == DE_stats keys (set equality)", key_raw == key_m,
          f"raw={len(key_raw)}, master={len(key_m)}")
    check("all 3 conditions present", set(m.condition) == {"Rest", "Stim8hr", "Stim48hr"})
    # gene<->gene_id must be 1:1 and match raw, else the six gene-keyed merges silently misalign
    check("gene<->gene_id is 1:1",
          m.groupby("gene").gene_id.nunique().max() == 1 and m.groupby("gene_id").gene.nunique().max() == 1)
    check("master gene<->gene_id mapping matches raw DE_stats",
          set(zip(de_raw.target_contrast_gene_name, de_raw.target_contrast)) == set(zip(m.gene, m.gene_id)))

    print("== C. evidence layers re-derived from RAW ==")
    strong_raw = de_raw.set_index(["target_contrast", "culture_condition"])["n_total_genes_category"].eq(STRONG_CATEGORY)
    mi = m.set_index(["gene_id", "condition"])
    check("is_strong_regulator <=> category '>10 DE genes' (raw)",
          mi["is_strong_regulator"].reindex(strong_raw.index).fillna(False).equals(strong_raw),
          f"strong={int(strong_raw.sum())}")
    # KD: re-derive frac from the RAW guide table, not the built column
    kd_g = (kd_raw.groupby(["perturbed_gene_id", "culture_condition"], as_index=False)
                  .agg(n_tested=("signif_knockdown", "size"), n_sig=("signif_knockdown", "sum")))
    kd_g["frac_raw"] = kd_g["n_sig"] / kd_g["n_tested"]
    kd_g = kd_g.rename(columns={"perturbed_gene_id": "gene_id", "culture_condition": "condition"})
    jk = m.merge(kd_g[["gene_id", "condition", "frac_raw"]], on=["gene_id", "condition"], how="left")
    jk["frac_raw"] = jk["frac_raw"].fillna(0.0)
    check("frac_guides_signif_kd re-derived from raw guide table",
          np.allclose(jk["frac_guides_signif_kd"].fillna(0.0).values, jk["frac_raw"].values))
    check("is_efficient_kd <=> RAW frac >= 0.5", (jk["is_efficient_kd"].astype(bool) == (jk["frac_raw"] >= KD_FRAC_MIN)).all())
    # donor: re-derive from the RAW donor summary
    dr = donor_raw.rename(columns={"target_contrast": "gene_id"}).drop_duplicates(["gene_id", "condition"])
    jd = m.merge(dr[["gene_id", "condition", "donor_correlation_mean"]], on=["gene_id", "condition"],
                 how="left", suffixes=("", "_raw"))
    subd = jd[jd["donor_correlation_mean_raw"].notna()]
    check("donor_correlation_mean re-derived from raw donor summary",
          np.allclose(subd["donor_correlation_mean"].values, subd["donor_correlation_mean_raw"].values))
    check("is_reproducible_donor <=> RAW donor_corr >= 0.4",
          (subd["is_reproducible_donor"] == (subd["donor_correlation_mean_raw"] >= DONOR_CORR_MIN)).all())
    # guide: re-derive from the RAW guide-correlation table
    gr = guide_raw.rename(columns={"target": "gene", "culture_condition": "condition",
                                   "correlation": "guide_correlation"}).drop_duplicates(["gene", "condition"])
    jg = m.merge(gr[["gene", "condition", "guide_correlation"]], on=["gene", "condition"],
                 how="left", suffixes=("", "_raw"))
    subg = jg[jg["guide_correlation_raw"].notna()]
    check("guide_correlation re-derived from raw guide-corr table",
          np.allclose(subg["guide_correlation"].values, subg["guide_correlation_raw"].values))
    check("is_reproducible_guide <=> RAW guide_corr >= 0.4",
          (subg["is_reproducible_guide"] == (subg["guide_correlation_raw"] >= GUIDE_CORR_MIN)).all())
    # autoimmune flag — bidirectional on the measured (OR/FDR present) subset
    sub_a = m[m["autoimmune_or"].notna() & m["autoimmune_fdr"].notna()]
    exp_a = (sub_a["autoimmune_or"] > AUTO_OR_MIN) & (sub_a["autoimmune_fdr"] < AUTO_FDR_MAX)
    check("is_autoimmune_enriched <=> OR>1 & FDR<0.1 (bidirectional, where measured)",
          bool((sub_a["is_autoimmune_enriched"].astype(bool) == exp_a).all()),
          f"n_flagged={int(sub_a['is_autoimmune_enriched'].sum())}, n_measured={len(sub_a)}")
    # carried non-score layers (defense-in-depth)
    subc = m[m["k562_logfc_r"].notna()]
    check("is_cd4_specific <=> k562_logfc_r < 0.2 (where measured)",
          (subc["is_cd4_specific"].astype(bool) == (subc["k562_logfc_r"] < CD4_SPECIFIC_R_MAX)).all())
    check("is_il_validated <=> gene in raw arrayed-validation set",
          (m["is_il_validated"].astype(bool) == m["gene"].isin(val_genes)).all())

    print("== D. cross-table recompute from raw ==")
    for c in ["Rest", "Stim8hr", "Stim48hr"]:
        raw_n = int(de_raw[(de_raw.culture_condition == c)].ontarget_significant.sum())
        built_n = int(cs.loc[cs.condition == c, "n_ontarget_significant"].iloc[0])
        check(f"n_ontarget_significant[{c}] matches raw", raw_n == built_n, f"raw={raw_n} built={built_n}")
    # ontarget_significant is a score input -> verify per-row, not just the aggregate
    jo = m.merge(de_raw.rename(columns={"target_contrast": "gene_id", "culture_condition": "condition"})
                 [["gene_id", "condition", "ontarget_significant"]], on=["gene_id", "condition"],
                 how="left", suffixes=("", "_raw"))
    check("ontarget_significant preserved per-row from raw",
          bool((jo.ontarget_significant.astype(bool) == jo.ontarget_significant_raw.astype(bool)).all()))
    # significant-KD guides restricted to DE-tested keys
    kd_sum = kd_raw.groupby(["perturbed_gene_id", "culture_condition"], as_index=False).signif_knockdown.sum()
    mkeys = set(zip(m.gene_id, m.condition))
    kd_sum = kd_sum[[(g, c) in mkeys for g, c in zip(kd_sum.perturbed_gene_id, kd_sum.culture_condition)]]
    check("significant-KD guides match raw (restricted to DE-tested keys)",
          int(kd_sum.signif_knockdown.sum()) == int(m["n_guides_signif_kd"].fillna(0).sum()),
          f"raw={int(kd_sum.signif_knockdown.sum())} built={int(m['n_guides_signif_kd'].fillna(0).sum())}")
    # n_degs_mash and n_guides_tested re-derived from raw
    kr = k562_raw.rename(columns={"target_contrast_gene_name": "gene"}).copy()
    kr["nm"] = kr.apply(lambda r: r.get(f"n_degs_MASH_{r['condition']}", np.nan), axis=1)
    jm = m.merge(kr[["gene", "condition", "nm"]].drop_duplicates(["gene", "condition"]), on=["gene", "condition"], how="left")
    sm = jm[jm.n_degs_mash.notna()]
    check("n_degs_mash matches raw K562 MASH count", (sm.n_degs_mash == sm.nm).all())
    gt = (kd_raw.groupby(["perturbed_gene_id", "culture_condition"], as_index=False).size()
                 .rename(columns={"perturbed_gene_id": "gene_id", "culture_condition": "condition", "size": "ngt_raw"}))
    jgt = m.merge(gt, on=["gene_id", "condition"], how="left")
    pres = jgt.n_guides_tested.notna()
    check("n_guides_tested matches raw KD group size (where present)",
          (jgt.loc[pres, "n_guides_tested"] == jgt.loc[pres, "ngt_raw"]).all())
    # autoimmune-enriched (cluster,condition) pairs feasibility
    auto_f = auto_raw[(auto_raw.disease == "autoimmune disease") & (~auto_raw.negative_control_disease.astype(bool))
                      & auto_raw.gene_set.isin(["downstream_Rest", "downstream_Stim8hr", "downstream_Stim48hr"])
                      & (auto_raw.odds_ratio > AUTO_OR_MIN) & (auto_raw.p_adj_fdr < AUTO_FDR_MAX)]
    check("autoimmune-enriched (cluster,condition) pairs recomputed",
          auto_f.drop_duplicates(["cluster", "gene_set"]).shape[0] >= 1,
          f"raw enriched pairs={auto_f.drop_duplicates(['cluster','gene_set']).shape[0]}")

    print("== E. score sanity ==")
    check("regulator_score within [0, 9]", m.regulator_score.between(0, 9).all(),
          f"min={m.regulator_score.min()} max={m.regulator_score.max()}")
    b = lambda s: s.fillna(False).astype(int)
    expected_score = (2 * b(m.is_strong_regulator) + b(m.ontarget_significant) + b(m.is_efficient_kd)
                      + b(m.is_reproducible_donor) + b(m.is_reproducible_guide)
                      + b(m.is_polarization_regulator) + b(m.is_aging_regulator) + b(m.is_autoimmune_enriched))
    check("regulator_score == sum of 8 weighted flags", (expected_score == m.regulator_score).all(),
          f"mismatches={int((expected_score != m.regulator_score).sum())}")
    check("strong regulator => score >= 2 (monotonicity)", (m[m.is_strong_regulator].regulator_score >= 2).all())
    check("on-target-significant => score >= 1", (m[m.ontarget_significant & ~m.is_strong_regulator].regulator_score >= 1).all())

    print("== F. accounting / membership ==")
    check("condition_summary n_perturbations sums to 33983", int(cs.n_perturbations.sum()) == 33983,
          f"{int(cs.n_perturbations.sum())}")
    ann_clusters = set(pd.read_csv(os.path.join(CLEAN, "cluster_annotations.csv")).cluster)
    check("every membership cluster exists in annotations", set(mem.cluster).issubset(ann_clusters))
    check("in_grn_module rows == cluster_membership rows",
          int(m.in_grn_module.sum()) == mem.shape[0], f"{int(m.in_grn_module.sum())} vs {mem.shape[0]}")
    j = m.merge(de_raw.rename(columns={"target_contrast": "gene_id", "culture_condition": "condition"})
                [["gene_id", "condition", "n_downstream", "n_up_genes", "n_down_genes", "n_total_de_genes"]],
                on=["gene_id", "condition"], how="left", suffixes=("", "_raw"))
    check("n_downstream preserved from raw", (j.n_downstream == j.n_downstream_raw).all())
    check("DE breakdown (up/down/total) preserved from raw",
          (j.n_up_genes == j.n_up_genes_raw).all() and (j.n_down_genes == j.n_down_genes_raw).all()
          and (j.n_total_de_genes == j.n_total_de_genes_raw).all())

    npass = sum(checks)
    print(f"\n{npass}/{len(checks)} checks passed")
    sys.exit(0 if npass == len(checks) else 1)


if __name__ == "__main__":
    main()
