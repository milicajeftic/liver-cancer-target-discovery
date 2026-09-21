# Finding Drug Targets in Liver Cancer from RNA-seq

> **Status:** In progress. 2 of 5 notebooks complete (quality control, differential expression).

**Question:** Which genes are overexpressed in hepatocellular carcinoma (HCC), linked to worse patient survival, *and* realistically druggable?

Hepatocellular carcinoma is the most common primary liver cancer and one of the leading causes of cancer death worldwide, with few effective drug options for advanced disease. This project uses public TCGA data to go from raw gene counts to a short, ranked list of candidate drug targets with a biological rationale for each.

**Full write-up:** [report.md](report.md)

## Pipeline

| Step | Notebook | What it does | Status |
|---|---|---|---|
| 1 | `01_data_download_and_qc.ipynb` | Download TCGA-LIHC counts + clinical data, sample QC, PCA, marker sanity checks | Done |
| 2 | `02_differential_expression.ipynb` | Tumor vs normal with PyDESeq2 (all samples + paired validation + sensitivity check), volcano plot, heatmap | Done |
| 3 | `03_pathway_enrichment.ipynb` | GSEA / over-representation (Hallmark, Reactome) with gseapy | Planned |
| 4 | `04_survival_analysis.ipynb` | Kaplan–Meier + Cox models for top upregulated genes | Planned |
| 5 | `05_druggability.ipynb` | Open Targets, DGIdb, Human Protein Atlas; final ranking | Planned |

## Data

- **Source:** [TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC) via the [UCSC Xena GDC hub](https://xenabrowser.net/)
- **Expression:** STAR gene counts (GENCODE v36)
- **Clinical:** overall survival, stage, grade
- **Samples after QC:** 371 primary tumors, 50 adjacent normal liver samples (50 patients with both)
- **Survival:** 365 patients with survival data, 131 deaths, median follow-up 19.6 months
- **Genes after filtering:** 22,107 of 60,660 (≥10 counts in ≥50 samples)

Data is not stored in the repo. Run `python -m src.data` to download it into `data/raw/`.

## Key results

### 1. Quality control

Normal liver samples cluster tightly, while tumors are highly heterogeneous, which is expected for HCC.

![PCA](figures/01_pca.png)

Known HCC markers behave as expected (GPC3 and AFP up, CYP2E1 and CYP1A2 down in tumors), confirming sample labels are correct.

![Markers](figures/01_marker_check.png)

**Outlier decisions:** 10 samples were flagged as >3 SD from their group centre on PCA. None were removed:
- 9 tumors: normal library sizes; they represent the extreme end of tumor heterogeneity, not technical failures.
- 1 normal (`TCGA-FV-A2QR-11A`): tumor markers within the normal range and mildly reduced CYP genes, consistent with diseased adjacent liver rather than tumor contamination. The removal rule was defined before inspecting the sample, and a sensitivity analysis confirmed it does not affect results (see below).

### 2. Differential expression

| Analysis | Up in tumor | Down in tumor |
|---|---|---|
| A · all samples (371 vs 50) | 4,496 | 1,307 |
| B · paired (50 patients) | 1,859 | 2,872 |
| **Robust (A and B)** | **1,842** | **1,261** |

_padj < 0.05 and |log2FC| > 1_

![Volcano plot](figures/02_volcano.png)

- The most significant upregulated genes are **cell-cycle / proliferation genes** (CDKN3, CDC25C, UBE2T, NUF2, CENPF, SKA1, TROAP).
- Known HCC genes were recovered without prior selection: **GPC3, TERT, IGF2BP1, HOXA13, REG3A, PLVAP, MAGEA1**.
- **Paired validation:** 99% of genes up in the paired analysis were also up in the main analysis (fold-change correlation r = 0.73).
- **Sensitivity check:** removing the borderline normal sample changed almost nothing (r = 0.9995, 196/200 top genes unchanged).

![Heatmap](figures/02_heatmap.png)

The 1,842 robust upregulated genes are carried forward to survival analysis.

_More results coming as the analysis progresses._

## Candidate targets

| Rank | Gene | log2FC | FDR | Survival HR | Druggability | Existing drugs | Rationale |
|---|---|---|---|---|---|---|---|
| _TODO_ | | | | | | | |

## Methods and design choices

- **Two complementary DE analyses:** all samples (maximum power, captures tumor heterogeneity) and paired patients only (controls for differences between individuals). Only genes significant in both are carried forward.
- **Thresholds set in advance:** padj < 0.05 and |log2FC| > 1; fold changes shrunk to reduce noise from low-count genes.
- **Pre-defined QC rules:** sample removal criteria decided before inspecting individual samples, with a sensitivity analysis for the borderline case.
- **Multiple testing:** Benjamini–Hochberg FDR throughout.
- **Survival (planned):** expression cut-offs chosen *before* looking at outcomes (median split), to avoid optimistic "best cut-off" p-values.

## Limitations

- Normal samples are *adjacent* non-tumor tissue, often from cirrhotic livers, so they are not truly healthy.
- Bulk RNA-seq mixes tumor, immune and stromal cells; some "tumor" signal may come from non-cancer cells.
- Some top genes by fold change are very high in only a subset of tumors; candidate ranking will also consider how many tumors overexpress each gene.
- mRNA expression ≠ protein level ≠ essentiality. Targets here are hypotheses that would need experimental validation.
- Median follow-up is 19.6 months; with 131 deaths survival analysis is reasonably powered, but long-term effects may be missed.

## Reproduce

```bash
# Option 1: conda
conda env create -f environment.yml
conda activate lihc-targets

# Option 2: venv + pip
python -m venv .venv
.venv\Scripts\activate        # Windows  (macOS/Linux: source .venv/bin/activate)
pip install -r requirements.txt

# then
python -m src.data            # downloads the data into data/raw/
jupyter lab                   # run notebooks in order
```

## Author

Milica Jeftic, Bioinformatics student, University of Primorska (UP FAMNIT)
[LinkedIn](#) · [GitHub](https://github.com/milicajeftic)