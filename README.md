# Finding Drug Targets in Liver Cancer from RNA-seq

>  Work in progress. Sections marked _TODO_ get filled in as the analysis is done.

**Question:** Which genes are overexpressed in hepatocellular carcinoma (HCC), linked to worse patient survival, *and* realistically druggable?

Hepatocellular carcinoma is the most common primary liver cancer and one of the leading causes of cancer death worldwide, with few effective drug options for advanced disease. This project uses public TCGA data to go from raw gene counts to a short, ranked list of candidate drug targets with a biological rationale for each.

![workflow](figures/00_workflow.png) <!-- TODO: add a simple workflow diagram -->

## Pipeline

| Step | Notebook | What it does | Status |
|---|---|---|---|
| 1 | `01_data_download_and_qc.ipynb` | Download TCGA-LIHC counts + clinical data, sample QC, PCA, marker sanity checks | ✅ |
| 2 | `02_differential_expression.ipynb` | Tumor vs normal with PyDESeq2 (paired design), volcano plot, heatmap | ⬜ |
| 3 | `03_pathway_enrichment.ipynb` | GSEA / over-representation (Hallmark, Reactome) with gseapy | ⬜ |
| 4 | `04_survival_analysis.ipynb` | Kaplan–Meier + Cox models for top upregulated genes | ⬜ |
| 5 | `05_druggability.ipynb` | Open Targets, DGIdb, Human Protein Atlas; final ranking | ⬜ |

## Data

## Data

- **Source:** [TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC) via the [UCSC Xena GDC hub](https://xenabrowser.net/)
- **Expression:** STAR gene counts (GENCODE v36)
- **Clinical:** overall survival, stage, grade
- **Samples after QC:** 371 primary tumors, 50 adjacent normal liver samples
- **Genes after filtering:** 22,107 of 60,660 (≥10 counts in ≥50 samples)

Data is not stored in the repo. Run `python -m src.data` to download it into `data/raw/`.

## Key results

## Key results

### 1. Quality control

Normal liver samples cluster tightly, while tumors are highly heterogeneous, which is expected for HCC.

![PCA](figures/01_pca.png)

Known HCC markers behave as expected (GPC3 and AFP up, CYP2E1 and CYP1A2 down in tumors), confirming sample labels are correct.

![Markers](figures/01_marker_check.png)

**Outlier decisions:** 10 samples were flagged as >3 SD from their group centre on PCA. None were removed:
- 9 tumors: normal library sizes; they represent the extreme end of tumor heterogeneity, not technical failures.
- 1 normal (`TCGA-FV-A2QR-11A`): tumor markers within the normal range and mildly reduced CYP genes, consistent with diseased adjacent liver rather than tumor contamination. The removal rule was defined before inspecting the sample. A sensitivity analysis is done in notebook 02.

_More results coming as the analysis progresses._

## Candidate targets

| Rank | Gene | log2FC | FDR | Survival HR | Druggability | Existing drugs | Rationale |
|---|---|---|---|---|---|---|---|
| _TODO_ | | | | | | | |

## Methods and design choices

- **Paired analysis:** where a patient has both tumor and normal samples, the patient is included in the model to control for individual variation.
- **Multiple testing:** Benjamini–Hochberg FDR throughout.
- **Survival:** expression cut-offs chosen *before* looking at outcomes (median split), to avoid optimistic "best cut-off" p-values.

## Limitations

- Normal samples are *adjacent* non-tumor tissue, often from cirrhotic livers, so they are not truly healthy.
- Bulk RNA-seq mixes tumor, immune and stromal cells; some "tumor" signal may come from non-cancer cells.
- mRNA expression ≠ protein level ≠ essentiality. Targets here are hypotheses that would need experimental validation.
- Median follow-up is only 19.6 months, which limits the statistical power of survival analysis.

## Reproduce

```bash
conda env create -f environment.yml
conda activate lihc-targets
python -m src.data          # downloads ~100 MB into data/raw/
jupyter lab                 # run notebooks in order
```

## Author

Milica Jeftic, Bioinformatics student, University of Primorska (UP FAMNIT)
[LinkedIn](#) · [GitHub](https://github.com/milicajeftic)
