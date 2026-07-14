#!/usr/bin/env python
"""
make_figures.py — figures for the CD4+ T cell master-regulator analysis.
Reads data/integrated/*.csv, writes figures/*.png (+ one graphml).
Run: <venv>/bin/python scripts/make_figures.py
"""
from __future__ import annotations
import os
import ast
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import networkx as nx

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root (portable)
INTEG = os.path.join(ROOT, "data", "integrated")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

# colorblind-safe palette
BLUE, AQUA, GREEN, YEL = "#0072B2", "#56B4E9", "#009E73", "#E69F00"
RED, VIOLET, MAGENTA, GREY = "#D55E00", "#7B3294", "#CC79A7", "#7F7F7F"
INK, GRID = "#222222", "#DDDDDD"
SEQ = plt.cm.Blues
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
CCOL = {"Rest": AQUA, "Stim8hr": BLUE, "Stim48hr": VIOLET}

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 150, "font.size": 10,
    "axes.edgecolor": INK, "axes.linewidth": 0.8, "axes.titlesize": 12,
    "axes.titleweight": "bold", "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.axisbelow": True, "figure.autolayout": False,
})


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print(f"  [fig] figures/{name}")


def load():
    m = pd.read_csv(os.path.join(INTEG, "regulator_master.csv"), low_memory=False)
    cs = pd.read_csv(os.path.join(INTEG, "condition_summary.csv"))
    cl = pd.read_csv(os.path.join(INTEG, "cluster_summary.csv"), low_memory=False)
    return m, cs, cl


# ---- fig1: trans-effect breadth (network hubs) ---------------------------
def fig_breadth(m):
    s = (m[(m.condition == "Stim8hr") & m.ontarget_significant]
         .nlargest(20, "n_downstream")[::-1])
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    vals = s.n_downstream.values
    ax.barh(s.gene, vals, color=SEQ(0.3 + 0.6 * vals / vals.max()), edgecolor=INK, linewidth=0.4)
    for y, v in enumerate(vals):
        ax.text(v + vals.max() * 0.01, y, f"{int(v):,}", va="center", fontsize=8)
    ax.set_xlabel("# downstream genes significantly affected (trans-effects, Stim8hr)")
    ax.set_title("Broadest CD4+ regulators — trans-effect hubs")
    ax.margins(x=0.12)
    save(fig, "fig1_regulator_breadth.png")


# ---- fig2: evidence-layer coverage + score distribution ------------------
def fig_coverage(m):
    layers = [("on-target KD sig.", "ontarget_significant"), ("efficient KD", "is_efficient_kd"),
              ("cross-donor repro.", "is_reproducible_donor"), ("cross-guide repro.", "is_reproducible_guide"),
              ("GRN module", "in_grn_module"), ("polarization reg.", "is_polarization_regulator"),
              ("aging reg.", "is_aging_regulator"), ("autoimmune-enriched", "is_autoimmune_enriched"),
              ("CD4-specific", "is_cd4_specific")]
    names = [n for n, _ in layers][::-1]
    counts = [int(m[c].fillna(False).sum()) for _, c in layers][::-1]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 5.2))
    a1.barh(names, counts, color=GREEN, edgecolor=INK, linewidth=0.4)
    for y, v in enumerate(counts):
        a1.text(v + max(counts) * 0.01, y, f"{v:,}", va="center", fontsize=8)
    a1.set_xlabel("# gene x condition rows flagged")
    a1.set_title("Evidence-layer coverage")
    a1.margins(x=0.14)
    sc = m.regulator_score.value_counts().sort_index()
    a2.bar(sc.index, sc.values, color=BLUE, edgecolor=INK, linewidth=0.4)
    for x, v in zip(sc.index, sc.values):
        a2.text(x, v + max(sc.values) * 0.01, f"{v:,}", ha="center", fontsize=7)
    a2.set_xlabel("regulator_score (0–9)")
    a2.set_ylabel("# gene x condition")
    a2.set_title("Convergent-evidence score distribution")
    a2.set_xticks(range(0, 10))
    save(fig, "fig2_evidence_coverage.png")


# ---- fig3: up/down DE per condition (diverging) --------------------------
def fig_direction(m):
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    up = [m.loc[m.condition == c, "n_up_genes"].sum() for c in CONDS]
    dn = [m.loc[m.condition == c, "n_down_genes"].sum() for c in CONDS]
    y = np.arange(len(CONDS))
    ax.barh(y, up, color=RED, edgecolor=INK, linewidth=0.4, label="upregulated")
    ax.barh(y, [-d for d in dn], color=BLUE, edgecolor=INK, linewidth=0.4, label="downregulated")
    for yi, (u, d) in enumerate(zip(up, dn)):
        ax.text(u, yi, f"  {u/1e6:.2f}M", va="center", fontsize=8)
        ax.text(-d, yi, f"{d/1e6:.2f}M  ", va="center", ha="right", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(CONDS)
    ax.axvline(0, color=INK, linewidth=0.9)
    ax.set_xlabel("total significant DE gene calls across all perturbations (down | up)")
    ax.set_title("Directional transcriptional burden by condition")
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.set_xticklabels([f"{abs(int(t))/1e6:.1f}M" for t in ax.get_xticks()])
    save(fig, "fig3_de_direction_by_condition.png")


# ---- fig4: top prioritized regulators (Stim8hr) --------------------------
def fig_top(m):
    s = m[m.condition == "Stim8hr"].nlargest(25, ["regulator_score"]).copy()
    s = s.sort_values(["regulator_score", "n_downstream"])[::1]
    fig, ax = plt.subplots(figsize=(8.2, 7.6))
    col = np.where(s.is_autoimmune_enriched.fillna(False), RED,
                   np.where(s.is_polarization_regulator.fillna(False), VIOLET, BLUE))
    ax.barh(range(len(s)), s.regulator_score, color=col, edgecolor=INK, linewidth=0.4)
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels(s.gene)
    for i, (_, r) in enumerate(s.iterrows()):
        tag = f"  ds={int(r.n_downstream):,}"
        ann = r.cluster_annotation if isinstance(r.cluster_annotation, str) else ""
        ax.text(r.regulator_score + 0.05, i, tag + (f" · {ann}" if ann and ann != "unknown" else ""),
                va="center", fontsize=7, color=INK)
    ax.set_xlabel("regulator_score (0–9)")
    ax.set_xlim(0, 9)
    ax.set_title("Top prioritized CD4+ regulators (Stim8hr)")
    ax.legend(handles=[Patch(color=RED, label="autoimmune-enriched module"),
                       Patch(color=VIOLET, label="polarization/activation reg."),
                       Patch(color=BLUE, label="other")],
              loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=3,
              frameon=False, fontsize=8)
    save(fig, "fig4_top_targets.png")


# ---- fig5: autoimmune-enriched GRN clusters ------------------------------
def fig_autoimmune(cl):
    s = cl[cl.best_autoimmune_or.notna() & (cl.min_autoimmune_fdr < 0.10)].copy()
    s = s.nlargest(15, "best_autoimmune_or").sort_values("best_autoimmune_or")
    lab = [f"{int(c)}: {a[:26] if isinstance(a, str) else '—'}"
           for c, a in zip(s.cluster, s.manual_annotation)]
    fig, ax = plt.subplots(figsize=(8.6, 6.0))
    ax.barh(range(len(s)), s.best_autoimmune_or, color=MAGENTA, edgecolor=INK, linewidth=0.4)
    ax.set_yticks(range(len(s)))
    ax.set_yticklabels(lab, fontsize=8)
    for i, (_, r) in enumerate(s.iterrows()):
        ax.text(r.best_autoimmune_or + 0.05, i, f"FDR={r.min_autoimmune_fdr:.1e}",
                va="center", fontsize=7)
    ax.set_xlabel("autoimmune-disease odds ratio (best per cluster)")
    ax.set_title("GRN modules enriched for autoimmune-disease genes")
    ax.margins(x=0.16)
    save(fig, "fig5_autoimmune_clusters.png")


# ---- fig6: regulator -> GRN module bipartite network ---------------------
def fig_network(m):
    s = m[(m.condition == "Stim8hr") & (m.regulator_score >= 8) & m.in_grn_module].copy()
    s = s[s.cluster_annotation.notna()]
    G = nx.Graph()
    for _, r in s.iterrows():
        reg, cl = r.gene, f"[{r.cluster_annotation}]"
        G.add_node(reg, kind="regulator", score=int(r.regulator_score))
        G.add_node(cl, kind="module")
        G.add_edge(reg, cl)
    if G.number_of_nodes() == 0:
        return
    nx.write_graphml(G, os.path.join(FIG, "regulator_module_network.graphml"))
    pos = nx.spring_layout(G, seed=7, k=0.9)
    fig, ax = plt.subplots(figsize=(12, 10))
    regs = [n for n, d in G.nodes(data=True) if d["kind"] == "regulator"]
    mods = [n for n, d in G.nodes(data=True) if d["kind"] == "module"]
    sc = [G.nodes[n]["score"] for n in regs]
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=GRID, width=0.8)
    nx.draw_networkx_nodes(G, pos, nodelist=mods, node_color=YEL, node_shape="s",
                           node_size=900, edgecolors=INK, linewidths=0.6, ax=ax)
    nc = nx.draw_networkx_nodes(G, pos, nodelist=regs, node_color=sc, cmap=plt.cm.viridis,
                                node_size=320, edgecolors=INK, linewidths=0.5, ax=ax)
    nx.draw_networkx_labels(G, pos, {n: n for n in mods}, font_size=8, font_weight="bold", ax=ax)
    nx.draw_networkx_labels(G, pos, {n: n for n in regs}, font_size=7, ax=ax)
    ax.set_title("High-confidence CD4+ regulators (score ≥ 8, Stim8hr) mapped to GRN modules")
    ax.axis("off")
    cb = fig.colorbar(nc, ax=ax, shrink=0.5)
    cb.set_label("regulator_score")
    save(fig, "fig6_regulator_module_network.png")


def main():
    m, cs, cl = load()
    fig_breadth(m)
    fig_coverage(m)
    fig_direction(m)
    fig_top(m)
    fig_autoimmune(cl)
    fig_network(m)
    print("[ok] figures written")


if __name__ == "__main__":
    main()
