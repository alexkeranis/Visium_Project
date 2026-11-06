# Methods

All parameters live in [`config.yaml`](../config.yaml). Code references point to
[`src/prostate_markers/`](../src/prostate_markers/).

## Data

Two 10x Genomics Visium FFPE datasets (Space Ranger 1.3.0, human probe set, 17,943 genes), one
donor each, tissue from Indivumed:

| Slide | 10x dataset | Spots under tissue |
|---|---|---|
| normal | Normal Human Prostate (FFPE) | 2,543 |
| cancer | Human Prostate Cancer, Adenocarcinoma with Invasive Carcinoma (FFPE), Gleason 7 | 4,371 |

Only the hires H&E image (about 2,000 px wide; one spot ≈ 14 px) is used. The full-resolution
image is not part of the bundled files.

## Quality control (`qc.py`)

Space Ranger's filtered matrix already contains only spots under tissue. We keep spots with
500-50,000 UMIs and at least 250 detected genes, which removes 1 normal and 2 cancer spots. The
FFPE probe set contains no mitochondrial genes, so mitochondrial fraction is not available as a
QC metric.

## Normalization (`normalize.py`)

Primary: counts per 10,000 followed by `log1p`, per spot (scanpy). Raw counts are kept in
`layers["counts"]`.

Sensitivity analyses:

- **Housekeeping scaling.** Counts divided by the per-spot mean of the configured housekeeping
  genes present in the probe set. GAPDH, RPLP0 and PPIA are absent, leaving ACTB and B2M.
  Both are lower in the cancer section and track depth within each slide
  (`results/tables/housekeeping_check.csv`). Housekeeping scaling is therefore not used as the
  primary method.
- **Epithelium-rich spots.** The two sections differ in composition (more luminal epithelium in
  the cancer section, more smooth-muscle stroma in the normal one). Fold changes are repeated
  using only the top third of spots in each slide by mean log-expression of EPCAM, KRT8 and
  KRT18.

## Gene effect sizes (`effects.py`)

For every gene, between the two slides:

- **log2 fold change** of pseudo-bulk expression: summed counts of the gene divided by all counts
  in the slide (per 10k), plus a pseudocount of 0.01, cancer over normal;
- **Cohen's d** on spot-level log1p(CP10k);
- **probability of superiority** P(X_cancer > X_normal), the Mann-Whitney U divided by
  n_cancer × n_normal (panel genes only);
- **detection fraction** per slide.

None of these is a test. With one section per class the spots are not independent replicates
of either condition, so we report no p-values for between-slide differences.

## Slide background (`background.py`)

Detected genes (≥10% of spots in either slide) are split into 20 quantile bins of mean
log-expression. For each panel gene we report the percentile of its |log2FC| (and |d|) among
all other genes in its bin. Every gene in the bin is used, so the result is deterministic. A gene
is *beyond background* at the 95th percentile or above.

For gene-set scores the reference is 500 random sets of detected genes of the same size (seeded).
We measure how well each set's mean expression separates the slides, as ROC AUC folded around 0.5.

## Literature direction (`panel.py`)

Each panel gene has an `expected` direction in cancer relative to benign tissue, and a PubMed ID
when a study supports it. `up`/`down` is assigned only when the cited study compares cancer
with benign tissue. `context` marks genes whose change depends on tumor subtype or grade.
`unknown` marks genes for which no single cancer-vs-benign reference is used here.

| Gene | Expected | Reference |
|---|---|---|
| AMACR | up | Rubin et al., JAMA 2002 · PMID 11926890 |
| FASN | up | Rossi et al., Mol Cancer Res 2003 · PMID 12939396 |
| MYC | up | Gurel et al., Mod Pathol 2008 · PMID 18567993 |
| NKX3-1 | down | Bowen et al., Cancer Res 2000 · PMID 11085535 |
| PTEN | down | Lotan et al., Clin Cancer Res 2011 · PMID 21878536 |
| TP63 | down (basal cells absent in cancer glands) | Weinstein et al., Mod Pathol 2002 · PMID 12481011 |
| KRT5 | down (basal cells absent in cancer glands) | Ali & Epstein, Am J Surg Pathol 2008 · PMID 18813120 |
| ERG, ETV1, ETV4 | context: up in ETS fusion-positive tumors | Tomlins et al., Science 2005 · PMID 16254181 |
| HOXB13 | context: overexpression prognostic within cancers | Zabalza et al., Oncotarget 2015 · PMID 25825985 |
| MKI67 | context: higher in aggressive disease | Fisher et al., Br J Cancer 2013 · PMID 23329234 |

Further background for the androgen-regulated kallikreins and for the panel as a whole:
Boyukozer et al., Croat Med J 2020 (PMID 33150763); Ligon et al., npj Precis Oncol 2026
(doi:10.1038/s41698-026-01368-3). Housekeeping-gene selection: Eisenberg & Levanon,
Trends Genet 2013 (PMID 23810203).

**Set scores.** `scanpy.tl.score_genes` (mean expression minus a control set matched on
expression level, seeded) over the genes expected *up* and the genes expected *down*. Scores are
computed on both slides together so they share one control set.

**Log-ratio.** KLK3 relative to AR, `log1p(KLK3) - log1p(AR)` on normalized data: an
androgen-response output relative to the receptor transcript.

## Tiers

| Criterion | Rule |
|---|---|
| beyond background | `bg_pct_log2fc` ≥ 95 |
| literature-concordant | expected `up`/`down` and observed in the same direction |
| direction robust | same sign under library-size, housekeeping and epithelium-only fold changes |

**A** = all three, **B** = two of three. Genes detected in fewer than 10% of spots in both slides
are not tiered. The tiers are a transparent ranking for follow-up, not a score with learned
weights.

## H&E image features (`texture.py`)

Square crop of side 2 × spot diameter (~27 px in the hires image), centered on each spot:

- hematoxylin and eosin optical density, the mean of the H and E channels after Ruifrok-Johnston
  color deconvolution (`skimage.color.rgb2hed`);
- mean and standard deviation of the grayscale crop;
- gray-level co-occurrence statistics (contrast, dissimilarity, homogeneity, energy,
  correlation) on 32 gray levels, symmetric and normalized, averaged over distances 1-2 px and
  angles 0°, 45°, 90°, 135°.

These features depend on staining, section thickness, scanner and quantization. They are not
invariant to these factors, so between-slide differences are reported only descriptively
(`image_between_slides.csv`).

## Within-slide concordance (`concordance.py`)

Within each slide separately: Spearman correlation between each molecular feature (panel genes
detected in ≥10% of spots in both slides, set scores, log-ratio) and each image feature.
A pair is *concordant* if the correlations share a sign and both have |rho| ≥ 0.2. Neighboring
spots are correlated, so the correlations are effect sizes; no p-values are computed.

## Reproducibility

- One seed (`config.yaml: seed`) seeds Python and NumPy. Each stochastic step (gene-set control
  genes, random gene sets) gets its own child seed from `numpy.random.SeedSequence`.
- Inputs are checked against `data/CHECKSUMS.sha256` before every run.
- Outputs have fixed names and fixed float formatting. Two clean runs give byte-identical tables
  and figures.
- `results/run_manifest.json` records the seed, config hash, git commit, package versions and
  input checksums of the last run.
