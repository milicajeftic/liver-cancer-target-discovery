# Identifying Candidate Drug Targets in Hepatocellular Carcinoma from TCGA RNA-seq Data

**Author:** Milica Jeftic, Bioinformatics, University of Primorska (UP FAMNIT)
**Status:** In progress. Section 1 (data and quality control) is complete.

---

## Contents

1. [Introduction](#1-introduction)
2. [Data](#2-data)
3. [Quality control](#3-quality-control)
4. Differential expression _(coming)_
5. Pathway enrichment _(coming)_
6. Survival analysis _(coming)_
7. Druggability assessment _(coming)_
8. [Limitations](#8-limitations)
9. [References](#9-references)

---

## 1. Introduction

Hepatocellular carcinoma (HCC) is the most common form of primary liver cancer and a leading cause of cancer-related death worldwide. Most patients are diagnosed at an advanced stage, when surgery is no longer an option, and the available systemic therapies extend survival only modestly. There is a clear need for new therapeutic targets.

This project uses publicly available RNA-seq data from The Cancer Genome Atlas (TCGA) to look for genes that meet three criteria:

1. **Overexpressed** in tumor compared with adjacent non-tumor liver
2. **Associated with worse overall survival** in patients
3. **Druggable**, meaning the protein belongs to a class that small molecules or antibodies can realistically target, or already has known compounds

Genes meeting all three are proposed as candidate drug targets. These are computational hypotheses and would require experimental validation.

---

## 2. Data

| Item | Details |
|---|---|
| Cohort | TCGA-LIHC (Liver Hepatocellular Carcinoma) |
| Access | UCSC Xena GDC hub |
| Expression | STAR gene-level read counts, GENCODE v36 annotation (60,660 genes) |
| Survival | Overall survival status (`OS`) and time in days (`OS.time`) |
| Clinical | 439 records, 83 variables, including pathologic stage, tumor grade, sex and age |

Xena provides counts as log2(count + 1). These were back-transformed to integer counts, which count-based methods such as DESeq2 require.

---

## 3. Quality control

_Notebook: `notebooks/01_data_download_and_qc.ipynb`_

### 3.1 Sample selection

Sample type was taken from the TCGA barcode (positions 14–15): `01` = primary tumor, `11` = solid tissue normal. Recurrent tumors and other sample types were excluded. Where a patient had more than one aliquot of the same sample type, one was kept (first vial alphabetically) so that no patient is counted twice.

| Group | Samples |
|---|---|
| Primary tumor | 371 |
| Adjacent normal liver | 50 |
| Patients with both tumor and normal | 50 (all normal samples have a matched tumor) |

### 3.2 Library size

Library size (total reads per sample) was checked to identify failed or unusual sequencing runs.

![Library sizes](figures/01_library_sizes.png)

| Group | Mean (M reads) | Median | Min | Max |
|---|---|---|---|---|
| Tumor | 48.3 | 47.9 | 18.6 | 119.7 |
| Normal | 54.0 | 52.8 | 25.0 | 76.7 |

No sample had a library size low enough to suggest a failed run. Two tumor samples had unusually high depth (~115–120 M reads). This is not a problem, because DESeq2 normalization corrects for sequencing depth.

### 3.3 Gene filtering

Genes with very low counts carry little information and increase the multiple-testing burden. A gene was kept if it had **at least 10 counts in at least 50 samples**, where 50 is the size of the smaller group (normals). This threshold keeps genes expressed only in normal tissue.

**Result:** 22,107 of 60,660 genes (36%) were kept. This proportion is typical, because many annotated genes (pseudogenes and tissue-specific genes) are not expressed in liver.

### 3.4 Principal component analysis

For visualization only, counts were normalized to log2 counts-per-million, and PCA was run on the 2,000 most variable genes.

![PCA](figures/01_pca.png)

- **PC1 explains 19.4%** of the variance and **PC2 explains 9.0%**.
- Normal samples form a tight cluster, consistent with the relatively uniform expression profile of non-tumor liver.
- Tumor samples are widely spread along PC1, reflecting the known molecular heterogeneity of HCC. Tumors closest to the normal cluster may be well-differentiated or have lower tumor purity. Tumors at the opposite extreme are likely the least liver-like.

### 3.5 Outlier review

Samples were flagged if their distance from their group's median position on PC1/PC2 was more than 3 standard deviations above the group's mean distance. Flagging was used to prompt manual review, not to remove samples automatically. The removal criterion was set **before** inspecting individual samples.

| Sample | Group | Library size (M) | Distance (z) | Decision |
|---|---|---|---|---|
| TCGA-BC-A10Q-01A | Tumor | 48.3 | 3.29 | Keep |
| TCGA-CC-5260-01A | Tumor | 32.2 | 4.04 | Keep |
| TCGA-CC-A3M9-01A | Tumor | 51.1 | 3.44 | Keep |
| TCGA-CC-A7II-01A | Tumor | 50.6 | 3.02 | Keep |
| TCGA-CC-A7IJ-01A | Tumor | 53.5 | 3.32 | Keep |
| TCGA-DD-A4NA-01A | Tumor | 53.5 | 3.16 | Keep |
| TCGA-ED-A7PX-01A | Tumor | 42.0 | 3.14 | Keep |
| TCGA-ED-A82E-01A | Tumor | 56.1 | 3.30 | Keep |
| TCGA-FV-A3I0-01A | Tumor | 46.0 | 3.47 | Keep |
| TCGA-FV-A2QR-11A | Normal | 50.6 | 3.23 | Keep |

**Flagged tumors (9).** All had normal library sizes, so there is no sign of technical failure. They sit at the extreme end of the tumor distribution along PC1 and most likely represent biological heterogeneity. Removing them would selectively exclude the most distinct tumors and could bias the tumor-vs-normal comparison.

**Flagged normal (TCGA-FV-A2QR-11A).** A normal sample that sits apart from the other normals could contain tumor tissue. Pre-defined criterion for removal: tumor markers (GPC3, AFP) above the maximum of all other normal samples **and** liver metabolism genes (CYP2E1, CYP1A2) reduced.

| Gene (log2 CPM) | This sample | Normal median | Normal max | Tumor median |
|---|---|---|---|---|
| GPC3 | 3.7 | 2.2 | 6.5 | 8.2 |
| AFP | 2.5 | 2.1 | 3.4 | 3.0 |
| CYP2E1 | 10.5 | 12.2 | 13.7 | 9.5 |
| CYP1A2 | 6.0 | 9.0 | 10.9 | 1.7 |

Tumor markers are within the normal range. CYP genes are mildly reduced but still within the range of normal samples. The removal criterion was not met, so the sample was kept. Its profile is consistent with diseased adjacent liver, since TCGA normal samples often come from cirrhotic or inflamed livers. A sensitivity analysis (differential expression with and without this sample) is planned in Section 4.

### 3.6 Marker gene validation

To confirm that sample labels and data processing are correct, four genes with well-established behavior in HCC were checked.

![Marker genes](figures/01_marker_check.png)

| Gene | Role | Expected in tumor | Observed |
|---|---|---|---|
| GPC3 | Established HCC biomarker | Higher | ✅ Strongly higher |
| AFP | Clinical serum marker, elevated in a subset of HCC | Higher in a subset | ✅ High in a subset, many tumor outliers |
| CYP2E1 | Hepatocyte metabolic enzyme | Lower | ✅ Lower, more variable |
| CYP1A2 | Hepatocyte metabolic enzyme | Lower | ✅ Strongly lower |

All four genes behave as expected, supporting the correctness of sample labeling.

### 3.7 Clinical and survival data

| Item | Value |
|---|---|
| Tumor samples with survival data | 365 |
| Deaths (events) | 131 |
| Median follow-up | 19.6 months |

Relevant clinical variables are available for later adjustment: `ajcc_pathologic_stage.diagnoses`, `tumor_grade.diagnoses`, `gender.demographic` and `age_at_index.demographic`.

### 3.8 Summary

The dataset passed quality control. No samples were removed. The final dataset for downstream analysis contains **371 tumor and 50 normal samples across 22,107 genes**. Marker genes confirm the sample labels, and PCA shows the expected separation between tumor and normal tissue alongside strong tumor heterogeneity.

---

## 8. Limitations

- **Adjacent normal is not healthy liver.** Normal samples come from non-tumor tissue of cancer patients, often with cirrhosis or hepatitis. Some tumor-vs-normal differences may be underestimated.
- **Bulk RNA-seq.** Each sample is a mixture of tumor, immune and stromal cells, so some signal may come from non-cancer cells, and tumor purity varies between samples.
- **Imbalanced groups.** There are 371 tumors vs 50 normals.
- **Short follow-up.** A median of 19.6 months limits the statistical power of survival analysis.
- **mRNA ≠ protein ≠ dependency.** High mRNA expression does not guarantee high protein levels or that the tumor depends on the gene.
- **Short follow-up.** Median follow-up is 19.6 months. With 131 deaths among 365 patients, survival analysis is still reasonably powered, but long-term effects may be missed.
---

## 9. References

1. Cancer Genome Atlas Research Network. Comprehensive and integrative genomic characterization of hepatocellular carcinoma. *Cell* 169, 1327–1341 (2017).
2. Goldman, M. J. et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology* 38, 675–678 (2020).
3. Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15, 550 (2014).
4. Muzellec, B., Teleńczuk, M., Cabeli, V. & Andreux, M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39, btad547 (2023).
