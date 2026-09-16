# Finding Drug Targets in Liver Cancer from RNA-seq

> 🚧 Work in progress. Sections marked _TODO_ get filled in as the analysis is done.

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

- **Source:** [TCGA-LIHC](https://portal.gdc.cancer.gov/projects/TCGA-LIHC) via the [UCSC Xena GDC hub](https://xenabrowser.net/)
- **Expression:** STAR gene counts (GENCODE v36)
- **Clinical:** overall survival, stage, grade
- **Samples:** primary tumors (`-01`) and solid tissue normal (`-11`); _TODO: fill in exact numbers after QC_

Data is not stored in the repo. Run `python -m src.data` to download it into `data/raw/`.

## Key results

_TODO_ — volcano plot, top pathways, Kaplan–Meier curves, final target table.

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
- _TODO: add what you find along the way_

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
