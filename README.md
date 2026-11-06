# Prostate Visium candidate markers

Exploratory, hypothesis-generating analysis of one normal and one cancer prostate section profiled
with 10x Genomics Visium (FFPE). We describe both tissues, then build a short, ranked list of
interpretable candidate markers from gene expression (literature-grounded) and H&E morphology.
Each candidate carries explicit evidence:

- its effect size;
- how its shift compares with the general slide-to-slide shift;
- whether it agrees with the published direction;
- whether its direction survives alternative normalizations;
- whether it tracks tissue morphology within each slide.

The candidates are meant as pre-specified inputs for testing on multi-patient cohorts.

**What this project does not do:** it trains no classifier, makes no diagnosis, and makes no claim
of generalization or patient-level significance. Spot labels mean *slide of origin*, not
pathologist-annotated tumor; the cancer section also contains benign glands.

## Main findings

| Gene | Tier | Change in cancer section | Published direction | Percentile among expression-matched genes | Epithelium-only log2FC |
|---|---|---|---|---|---|
| **AMACR** | A | up, log2FC 5.8 | up (PMID 11926890) | 99.6 | 6.0 |
| **TP63** | B | down, −1.4 | down, basal-cell loss (PMID 12481011) | 93.4 | −3.1 |
| **KRT5** | B | down, −1.3 | down, basal-cell loss (PMID 18813120) | 91.3 | −3.1 |
| **PTEN** | B | down, −1.5 | down (PMID 21878536) | 93.7 | −2.1 |
| **MYC** | B | up, 0.9 | up (PMID 18567993) | 71.1 | 1.1 |
| **ERG** | B | up, 3.6 | fusion-positive tumors only (PMID 16254181) | 99.4 | 4.3 |

FKBP5, KLK2 and HOXB13 also reach tier B (beyond background and robust) but have no
cancer-vs-benign reference in the panel. Full table:
[`results/tables/candidate_markers.csv`](results/tables/candidate_markers.csv).
Tier rules: [docs/methods.md](docs/methods.md#tiers).

- **AMACR** is the strongest candidate. It meets all three criteria, and its rise exceeds 99.6% of
  genes with similar expression.
- **The basal-cell genes TP63 and KRT5** fall in the published direction, and the drop doubles
  within epithelium-rich spots. That points to loss of the basal layer rather than just less
  stroma.
- **The "up" score** (AMACR, FASN, MYC) separates the slides better than 99.8% of random 3-gene
  sets.
- **The sections differ in composition.** The cancer section is richer in luminal epithelium,
  which moves every luminal gene up. NKX3-1, expected to be lost in cancer, rises with them; the
  increase shrinks when only epithelium-rich spots are compared.
- **Luminal genes follow gland morphology in both sections.** Luminal-epithelial genes (NKX3-1,
  HOXB13, TMPRSS2, KLK2/3, CDH1) rise with hematoxylin density and gray-level spread (gland
  lumens) and fall with eosin density. This is a tissue-architecture link, independent of slide
  identity, and not a cancer signal.
- **AMACR follows morphology in the cancer section only.** There it tracks hematoxylin density
  (rho 0.55) and runs against eosin (rho −0.55), marking nucleus-dense, stroma-poor regions. It
  does not do so in the normal section. A cancer-only link cannot be separated from slide identity
  here, so it is a lead to test on annotated sections.

<p align="center">
  <img src="results/figures/02_background.png" width="49%" alt="Panel genes against all detected genes">
  <img src="results/figures/02_panel_directions.png" width="40%" alt="Fold change under three normalizations">
</p>
<p align="center">
  <img src="results/figures/02_AMACR_spatial.png" width="90%" alt="AMACR expression over both sections">
</p>

## Approach

```
raw Visium (h5 + hires H&E)
  └─ prepare ─ checksum inputs · spot QC · CP10k + log1p
       └─ image ─ per-spot H&E stain density + co-occurrence features
            └─ markers ─ full-transcriptome effect sizes
                 ├─ slide background (expression-matched genes, random gene sets)
                 ├─ literature direction (config panel with PMIDs)
                 ├─ direction under housekeeping scaling and epithelium-only spots
                 ├─ within-slide gene ↔ image correlations
                 └─ candidate_markers.csv (tiers A / B)
                      └─ report ─ notebooks 01-03 → results/figures
```

**Why no classifier?** With one section per class, any classifier would learn to tell the two
slides apart (donor, fixation, staining, sequencing depth). Its accuracy would say nothing about
cancer. We instead ask whether specific, named genes shift more than the slide background,
in the direction the literature predicts, and robustly to normalization choices.

The notebooks walk through each step:

1. [`01_data_qc`](notebooks/01_data_qc.ipynb): QC, sequencing depth, tissue composition.
2. [`02_gene_markers`](notebooks/02_gene_markers.ipynb): background, housekeeping check,
   candidate table, set scores, spatial maps.
3. [`03_morphology_concordance`](notebooks/03_morphology_concordance.ipynb): H&E features and
   within-slide concordance.

Details and references: [docs/methods.md](docs/methods.md).

## Quickstart

Python 3.12. From the repository root:

```bash
python -m pip install -r requirements.txt
python -m pip install -e . --no-deps
python scripts/run_pipeline.py            # all stages, about 2 minutes
python scripts/run_pipeline.py --stage markers   # or a single stage
pytest -q
```

Stages: `prepare`, `image`, `markers`, `report` (`report` executes the notebooks).
`--skip-report` runs everything except the notebooks.

The run writes:

| Path | Tracked | Contents |
|---|---|---|
| `results/tables/` | yes | candidate table, QC summary, housekeeping check, image concordance, random-set reference |
| `results/figures/` | yes | figures saved by the notebooks |
| `results/run_manifest.json` | yes | seed, config hash, git commit, package versions, input checksums |
| `data/processed/` | no | normalized `.h5ad` per slide, per-spot features, all-gene effects |

All parameters, including the seed and the gene panel, are in [`config.yaml`](config.yaml).
Two clean runs produce byte-identical tables and figures.

## Scope

- **One section per class, from two donors.** Every between-slide difference mixes cancer
  biology with donor, tissue handling, staining and sequencing depth. Spots are not independent
  replicates, so no p-values are reported for between-slide differences.
- **Spot labels are slide labels.** Visium spots (55 µm) mix several cells. The cancer section
  contains benign glands and stroma, so *cancer* means "from the cancer section".
- **Composition.** The sections differ in epithelium/stroma content. Epithelium-only fold changes
  reduce, but do not remove, this effect.
- **Housekeeping genes.** Only ACTB and B2M of the configured set are in the probe set, and both
  shift between slides. Housekeeping scaling is a sensitivity check only.
- **Image resolution.** Image features use the hires image (about 14 px per spot). Full-resolution
  images would support finer morphology.

## Proposed validation on larger cohorts

1. Pre-specify the tier A/B genes and the up/down set scores from this repository, and do not
   reselect features on the new data.
2. Use many patients per class, with splits grouped by patient (and by slide or batch).
3. Use pathologist-annotated tumor and benign regions, so that spot-level labels mean tumor.
4. Apply stain normalization before any between-slide comparison of image features.
5. Test ERG together with a fusion assay (TMPRSS2-ERG), and NKX3-1 within annotated tumor glands.

## Data and attribution

Both datasets are from 10x Genomics and are distributed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). They are included in `data/` with
SHA-256 checksums in [`data/CHECKSUMS.sha256`](data/CHECKSUMS.sha256).

- *Normal Human Prostate (FFPE)*, Space Ranger 1.3.0, 10x Genomics (2021).
  <https://www.10xgenomics.com/datasets/normal-human-prostate-ffpe-1-standard-1-3-0>
- *Human Prostate Cancer, Adenocarcinoma with Invasive Carcinoma (FFPE)*, Space Ranger 1.3.0,
  10x Genomics (2021).
  <https://www.10xgenomics.com/datasets/human-prostate-cancer-adenocarcinoma-with-invasive-carcinoma-ffpe-1-standard-1-3-0>

## Repository layout

```
config.yaml              parameters, seed, gene panel with expected directions and PMIDs
data/                    raw 10x inputs + checksums (processed/ is generated, not tracked)
src/prostate_markers/    io, qc, normalize, effects, background, panel, texture, concordance, viz
scripts/run_pipeline.py  single entry point
notebooks/               01-03, executed by the report stage
results/                 tables, figures, run manifest
docs/methods.md          methods and references
tests/                   unit tests for the statistical helpers
```

## License

Code: [MIT](LICENSE). Data: CC BY 4.0, 10x Genomics. To cite this repository, see
[`CITATION.cff`](CITATION.cff).
