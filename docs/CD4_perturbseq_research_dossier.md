# Research dossier — CD4⁺ T cell Perturb-seq hackathon

**Purpose.** Arm the hackathon deliverable for the *Primary Human CD4⁺ T Cell Perturb-seq*
dataset (Zhu, Dann, … Pritchard, Marson; bioRxiv, 24 Dec 2025) with the relevant published
research, mapped onto the integrated-master-table deliverable style you already proved out in
`analysis.md` (the SARS-CoV-2 host-factor drug-repurposing analysis). It is built from a
5-angle deep-research sweep (23 sources fetched, 25 claims put through 3-vote adversarial
verification; 23 confirmed, 2 refuted). Confidence tags and the two **do-not-cite** refutations
are carried through honestly below.

---

## 0. TL;DR — what the best deliverable is, and the research that enables it

Your winning move is the **same play as the SARS-CoV-2 analysis, in a new domain**: integrate
the dataset's already-shipped functional-genomics layers into **one scored master table of
CD4⁺ T cell "master regulators,"** rank them, and annotate each with reproducibility, biology
(polarization / cytokine / aging), disease-genetics support, and druggability — with honest
scoring and caveats.

The reason this is tractable in a hackathon window: **the authors ship the hard parts
pre-computed** (genome-wide DESeq2 per-perturbation stats, pseudobulk, HDBSCAN trans-effect
modules, Th1/Th2 validation, aging signature, autoimmune enrichment, K562 comparison) and the
**processed data is public with no credentials** (`aws s3 ls --no-sign-request
s3://genome-scale-tcell-perturb-seq/marson2025_data/`). You do not need to touch the 22M raw
cells. Your value-add is the **integration + prioritization + honest scoring**, exactly as
before.

Three references do most of the work:
1. **The dataset's own paper** (angle 1) — tells you the intended headline biology and gives
   you the abstract-level claims to reproduce/extend from the shipped tables.
2. **Freimer et al. Nat Genet 2022** (angle 4/5) — the closest prior template: build a
   disease-enriched T cell gene-regulatory network from perturbation **trans-effects**, then
   link it to GWAS risk under an omnigenic model. This is your deliverable's backbone method,
   already validated in the exact cell type.
3. **Schmidt et al. Science 2022** (angle 4) — the known-vs-novel benchmark: which regulators a
   genome-scale primary-T-cell screen *should* recover, so you can label your hits
   "known / novel" credibly.

---

## 1. The flagship dataset & its paper (angle 1) — the primary source

**Zhu R, Dann E, Yan J, Reyes Retana J, Goto R, Guitche RC, Petersen LK, Ota M, Pritchard JK,
Marson A.** *"Genome-scale perturb-seq in primary human CD4⁺ T cells maps context-specific
regulators of T cell programs and human immune traits."* **bioRxiv, 24 Dec 2025.** DOI
`10.64898/2025.12.23.696273`. [confidence: **high**, 3-0]

What it establishes that your deliverable should anchor on:

- **Design/scale.** Genome-scale CRISPRi, probe-based Perturb-seq: knock down *all expressed
  protein-coding genes* across **~22M primary human CD4⁺ T cells, 4 donors, 3 conditions
  (Rest, Stim8hr, Stim48hr).** [3-0]
- **Central thesis = context-specificity.** "Active regulators and the gene programs they
  control change dramatically across stimulation conditions." **Design implication: your master
  table must be stratified by Rest / Stim8hr / Stim48hr, not pooled.** [3-0]
- **Headline biology outputs** (the things to surface): nominates **regulators of Th1 and Th2
  polarization** and of **age-related (aging) CD4⁺ phenotypes**; identifies **novel regulators
  of cytokine production**; **maps gene-regulatory networks**; and **models T cell states seen
  in population-scale atlases**. It also **implicates context-specific regulatory pathways in
  autoimmune-disease risk.** [3-0]
- **Companion model — scLDM.CD4.** Transformer-autoencoder + conditional flow-matching
  Diffusion Transformer, trained on **~14.5M cells / 3,699 HVGs** with donor/timepoint/guide
  metadata. Intended use: **synthesize single-cell profiles under single-gene knockdowns** and
  **rank candidate perturbations toward a desired transcriptomic effect.** Model cards: CZI
  Virtual Cell Models, HuggingFace `biohub/scldm_cd4`, arXiv `2511.02986`. [3-0]
  - ⚠️ **Honest caveat (important for scoring):** these are **developer intended-use
    statements, not independently benchmarked accuracy.** The cards restrict it to single-gene
    knockdowns and to conditions present in training (the same 4 donors / 3 timepoints),
    prohibit generating unseen gene combinations, and state outputs "require review and
    validation." **Cite it as a prioritization aid, not a proven predictor** — and if you use
    its rankings, benchmark against a trivial baseline (pseudobulk mean-shift / linear model).

**Ready-to-use products & pipeline** (GitHub `emdann/GWT_perturbseq_analysis_2025`) [3-0]:
processed cell-level counts, pseudobulk counts, and **DESeq2 DE estimates** on public S3
(`--no-sign-request`). The official pipeline is **8 sequential modules whose names map onto
your angles**: `3_DE_analysis`, `4_polarization_signatures`, `5_cytokine_regulators`,
`6_functional_interaction`, `7_1k1k_analysis` (OneK1K population eQTL/GWAS cohort),
`8_lymphocyte_counts_LoF` (heritable blood/immune traits). **Use these as reference code so you
don't re-derive DE from 22M cells.**

---

## 2. Recommended deliverable design — the CD4⁺ "master-regulator" scored table

Direct analog of your `host_factor_master.csv`. One row per **perturbed regulator × culture
condition** (the DE_stats obs schema is already this shape: `n_obs = 33,983`). Suggested scored
columns, each backed by a **shipped file**, so the whole thing is reproducible in-window:

| Layer (→ analog of SARS-CoV-2 table) | Column(s) | Source file in the dataset |
|---|---|---|
| **Regulatory breadth** (≈ "interaction burden") | `n_downstream`, `n_total_de_genes` | `GWCD4i.DE_stats.h5ad` / `DE_stats.suppl_table.csv` |
| **On-target validity** (QC gate) | `ontarget_significant`, `low_target_gex`, `signif_knockdown` | DE_stats + `guide_kd_efficiency.suppl_table.csv` |
| **Reproducibility** (≈ your "honesty" axis) | `guide_correlation_signif`, `donor_correlation_hits_mean` | DE_stats (cross-guide, cross-donor) |
| **Polarization role** | Th1/Th2 z-scores + **arrayed CRISPRi/flow validation** | `Th2_Th1_polarization…csv`, `Th1Th2_validation_summary…csv`, `polarization_prediction…coefficients.csv` |
| **Aging role** | aging-signature coefficients | `CD4T_aging_signature_DE…csv`, `aging_prediction…coefficients.csv` |
| **GRN module membership** | HDBSCAN cluster + `sign_coherence` | `clustering_downstream_genes.csv` |
| **Disease genetics** (≈ "pan-viral" reuse axis) | autoimmune odds-ratio / FDR by cluster | `cluster_autoimmune_enrichment_results.csv` |
| **Cell-type specificity** | K562-vs-CD4 logFC correlation | `K562_comparison.suppl_table.csv` |
| **Druggability** (the repurposing headline) | existing modulators / tractability | **external** — DGIdb / Open Targets / ChEMBL (see §5) |

**Scoring heuristic (keep it transparent, like your `evidence_score`):** gate on on-target KD
significance + reproducibility (cross-guide & cross-donor), then rank by regulatory breadth,
and **award bonus points for convergence** — a regulator that is *also* a validated polarization
hit *and* sits in an autoimmune-enriched cluster *and* has a drug is your "HDAC2/LARP1/PRKACA"
equivalent. As in `analysis.md`, **report the score as a heuristic, not a model**, and include
an honesty note wherever a genetically-important regulator scores low only because a layer
(e.g., druggability) is missing — the TOMM70 move.

The two most defensible "novel angles" for a hackathon: (a) **cross-condition regulators** —
genes whose effect flips or emerges only on stimulation (the paper's own thesis, quantifiable
via MASH); (b) **disease-anchored regulators** — cluster members that are autoimmune-enriched
*and* druggable but *not* in the prior-screen benchmark set (i.e., genuinely new therapeutic
leads).

---

## 3. Methods playbook (angle 2)

- **MASH — multivariate adaptive shrinkage.** Urbut, Wang, Carbonetto & Stephens, **Nature
  Genetics 2019** (`s41588-018-0268-8`; bioRxiv `096552`). Empirical-Bayes method that *learns*
  the correlation/sharing structure among conditions and exploits it to raise power, improve
  effect estimates, and **quantify effect-size heterogeneity** — rather than a binary
  shared-vs-specific call. **This is the right tool to share DESeq2 log2FCs across
  Rest/Stim8hr/Stim48hr and to define "context-specific" regulators quantitatively.** [3-0]
- **CRISPRi knockdown-efficiency QC.** Replogle et al., **eLife 2022** (`81856`): **dual-sgRNA
  median KD 82% vs single-sgRNA Dolcetto 65%** (Mann-Whitney p=2.4e-7). Use as the QC
  expectation when reading per-perturbation effect sizes / the shipped `guide_kd_efficiency`
  table. [3-0] ⚠️ **Do NOT cite** the related claim that "Zim3-dCas9 is unambiguously the
  optimal effector with no downsides" — **refuted 1-2.**
- **Genome-scale Perturb-seq analysis playbook.** Replogle, Weissman et al., **Cell 2022**,
  185:2559–2575 (`10.1016/j.cell.2022.05.013`; resource: gwps.wi.mit.edu). Establishes the
  standard moves your deliverable can borrow: **minimal-distortion embedding, clustering
  perturbations into pathways/modules, and trans-effect networks.** [reference is standard; note
  the live gwps.wi.mit.edu resource page returned no machine-extractable text during the sweep,
  so cite the paper, not the site.]
- **Causal GRN inference from perturbation data — LLCB (Linear Latent Causal Bayes).** *Cell
  Genomics 2024* (bioRxiv 2023-09-17). Bayesian causal structure-learning that infers directed
  GRNs directly from CRISPR perturbation in **primary human CD4⁺ T cells** (84 genes: 30 IEI
  TFs, 30 background TFs, 24 IL2RA regulators → **211 directed edges**). A directly citable
  blueprint if you want to go beyond correlation modules to **directed** regulator→target
  edges. [surfaced in fetch; not in the 3-0 verified core — confirm authorship/details before
  citing.]
- **Network master-regulator inference — ARACNe → VIPER.** *Cancer Cell 2023*
  (`S1535-6108(23)00129-0`) applied ARACNe network reconstruction + VIPER protein-activity
  inference to nominate **17 master regulators** of a tumor-infiltrating Treg state. A
  well-established alternative "master regulator" framing you can name-check or use as a
  cross-check on the perturbation-derived hits. [surfaced in fetch; not in verified core.]

---

## 4. Known-regulator benchmark set (angles 3 + 4) — to label hits "known vs novel"

**Canonical CD4⁺ T-helper master TFs the screen should recover as positive controls** (from
standard immunology references — O'Shea/Paul lineage reviews, *Nat Rev Immunol* 2013 `nri3321`;
Zhu & Paul; STAT reviews): **Th1 → T-bet/TBX21 (STAT1/STAT4, IFNG); Th2 → GATA3 (STAT6,
IL4/IL5/IL13); Th17 → RORγt/RORC (STAT3, IL17); Treg → FOXP3 (IL2/IL2RA); Tfh → BCL6.** [these
are textbook; angle-3 coverage in the verified set only nailed TBX21 & GATA3 via the Schmidt
benchmark, so **supply the RORC/FOXP3/BCL6/STAT panel from a standard review** rather than from
the sweep's claims.]

**Prior CRISPR functional-genomics screens in primary human T cells** — your "already known"
list, so novel hits stand out:

- **Shifrut et al., Cell 2018** (`S0092-8674(18)31333-3`; PMID 30449619). **SLICE** (sgRNA
  lentiviral infection + Cas9 protein electroporation) — genome-wide LOF proliferation screen,
  **77,441 sgRNAs / 19,114 genes.** The foundational primary-T-cell screen. [3-0] ⚠️ **Do NOT
  cite the specific negative-regulator list** "DGKZ/RhoH/SOCS1/CBLB/TCEB2/RASA2" as-is —
  **refuted 1-2**; verify individual hit names against the paper before naming them.
- **Schmidt, Steinhart et al., Science 2022** (`abj4008`; PMC9307090). Genome-wide **paired
  CRISPRa + CRISPRi** screens in primary human CD4⁺/CD8⁺ T cells (**>18,800 genes, >112,000
  sgRNAs**) for stimulation-induced **IL-2 and IFN-γ**. Recovered canonical regulators —
  positive: **VAV1, CD28, LCP2/SLP-76, TBX21**; negative: **MAP4K1, CBLB, GATA3** — plus novel
  **FOXQ1, APOBEC3A/D/C, EMP1**. Screen hits showed **enriched heritability for immune/
  autoimmune traits** vs non-immune (stratified LDSC). **This is your best known-vs-novel
  benchmark and it demos the GWAS-heritability-enrichment method you'll reuse.** [3-0]
- **Freimer et al., Nat Genet 2022** (`s41588-022-01106-y`; PMC10035359). SLICE + FACS for
  regulators of **IL2RA / IL-2 / CTLA4** (51/66/59 hits; 117 unique, 10 shared). **24 top IL2RA
  regulators form a dense feedback GRN** (each KO alters a median 9.5 of the other 24), the
  network is **significantly enriched for Mendelian + GWAS immune-disease genes**, and **17
  MS-associated SNPs fall in KO-altered ATAC peaks** — framed by the **omnigenic model.**
  **This is the single closest template for your trans-effect-GRN → disease-gene deliverable.**
  [3-0]

---

## 5. Perturbation → disease → druggability (angle 5) — the therapeutic payoff

- **In-dataset disease linkage (use first):** `cluster_autoimmune_enrichment_results.csv`
  (Fisher-exact autoimmune-gene enrichment per HDBSCAN cluster, with negative-control diseases)
  and pipeline modules `7_1k1k_analysis` (**OneK1K** population-scale immune-cell eQTL/GWAS
  cohort) and `8_lymphocyte_counts_LoF`. These let you connect clusters → disease directly from
  shipped outputs. [3-0]
- **Perturb-seq + eQTL + GWAS integration — "Mr. PEG."** *medRxiv, 5 Jan 2026*
  (`10.64898/2026.01.05.26343421`). Statistical framework jointly integrating perturbational
  (Perturb-seq) screens, eQTL, and GWAS summary stats to identify **trans-acting mediating
  genes** for complex traits — **the exact analytic bridge** from CD4⁺ trans-effects to
  disease-gene prioritization. Likely a companion to the flagship paper; **confirm before
  citing** (surfaced in fetch, not in verified core).
- **Treg QTL colocalization template.** Bossini-Castillo et al., *Cell Genomics 2022*
  (PMC9010307): eQTL + chromatin-QTL in **Treg from 124 individuals**, **133 colocalizing loci**
  where immune-disease GWAS variants overlap Treg regulatory activity. A concrete GWAS-to-gene
  mapping precedent in CD4⁺ regulatory T cells. [fetch-level.]
- **Causal-gene / druggability resource.** Open Targets Genetics — Mountjoy et al., **Nat Genet
  2021** (`s41588-021-00945-5`): fine-maps and prioritizes causal genes across **133,441 GWAS
  loci** with colocalization across 92 cell types. Pair with **DGIdb / ChEMBL** to annotate
  which top regulators have **approved drugs or clinical-stage modulators** — this is the
  column the verified set does *not* cover, so **build druggability from these external
  resources** (as `analysis.md` did with its 69-compound table). [fetch-level.]

---

## 6. Honest caveats (carry these into the deliverable)

- **Preprint, not peer-reviewed.** The flagship paper is a Dec-2025 bioRxiv preprint; "novel"
  is the authors' own framing. During the sweep, direct WebFetch of the bioRxiv page returned
  HTTP 403, so several confirmations rest on multi-host mirrors (CZI, SSRN, GitHub) rather than
  the rendered page — the facts are consistent across ≥2 primary hosts, but treat exact
  in-text numbers as "confirm in full text."
- **Two refuted claims — do not cite:** (a) Zim3-dCas9 as unambiguously optimal effector (1-2);
  (b) the specific Shifrut negative-regulator list DGKZ/RhoH/SOCS1/CBLB/TCEB2/RASA2 (1-2).
- **scLDM.CD4 is unbenchmarked** (see §1 caveat) — prioritization aid, not proven predictor.
- **Coverage gaps the sweep could not fill from the abstract alone** (need full-text/supplement
  or the actual data): the paper's **specific named novel cytokine/polarization/aging
  regulators**; the concrete **HDBSCAN module → autoimmune-disease** mappings; the **K562-vs-CD4**
  specifics; and per-regulator **druggability**. These are exactly where your hackathon
  value-add lives — they're computable from the shipped tables.

---

## 7. Open questions to resolve from the full text / the data

1. Which exact TFs/cytokine controllers does the preprint name as its **novel** Th1/Th2, aging,
   and cytokine-production regulators? (populate the master table's "novel" flag)
2. What are the concrete **HDBSCAN downstream-gene modules**, and which **autoimmune diseases**
   are enriched in which clusters? (`clustering_downstream_genes.csv` +
   `cluster_autoimmune_enrichment_results.csv`)
3. Which top regulators already have **approved drugs / clinical modulators**? (DGIdb / Open
   Targets / ChEMBL)
4. How does **scLDM.CD4** perform vs simple baselines on held-out data? (benchmark before trust)
5. What does **K562-vs-CD4** reveal about cell-type-specific vs shared regulators?
   (`K562_comparison.suppl_table.csv`)

---

## 8. Citation quick-list (verified core in **bold**)

Primary source & resources:
- **Zhu, Dann, … Pritchard, Marson. bioRxiv 2025.** DOI 10.64898/2025.12.23.696273 — flagship.
- **scLDM.CD4** — CZI Virtual Cell Models model card; HF `biohub/scldm_cd4`; arXiv 2511.02986.
- **GitHub `emdann/GWT_perturbseq_analysis_2025`** — pipeline + public S3 data.

Methods:
- **Urbut, Wang, Carbonetto, Stephens. MASH. Nat Genet 2019** (s41588-018-0268-8).
- **Replogle et al. eLife 2022** (81856) — CRISPRi KD-efficiency QC.
- Replogle, Weissman et al. Cell 2022, 185:2559–2575 — genome-scale Perturb-seq playbook.
- LLCB — Cell Genomics 2024 (causal GRN from CD4⁺ perturbation). *confirm details.*
- ARACNe/VIPER master-regulator analysis — Cancer Cell 2023 (S1535-6108(23)00129-0).

Benchmark screens & immunology:
- **Shifrut et al. Cell 2018** (PMID 30449619) — SLICE genome-wide KO.
- **Schmidt, Steinhart et al. Science 2022** (abj4008) — CRISPRa/i cytokine screens.
- **Freimer et al. Nat Genet 2022** (s41588-022-01106-y) — IL2RA GRN → disease (template).
- O'Shea/Paul lineage-TF reviews — Nat Rev Immunol 2013 (nri3321) and related.

Disease / druggability:
- Mr. PEG — medRxiv 2026 (10.64898/2026.01.05.26343421) — Perturb-seq+eQTL+GWAS. *confirm.*
- Bossini-Castillo et al. Cell Genomics 2022 (PMC9010307) — Treg QTL colocalization.
- Mountjoy et al. Open Targets Genetics, Nat Genet 2021 (s41588-021-00945-5); DGIdb; ChEMBL.

---
*Generated by a 5-angle deep-research sweep (23 sources, 25 claims adversarially verified,
23 confirmed / 2 refuted). "Verified core" = passed 3-0 verification against primary sources;
"fetch-level" items were extracted from primary sources but not part of the verified core —
confirm before citing exact figures.*
