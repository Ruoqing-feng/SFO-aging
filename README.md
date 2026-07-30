# SFO-aging

Analysis code associated with the manuscript:

**“Blood-accessible borders govern white matter ageing”**  
Liu *et al.* (2026)

## Scope of this repository

No novel standalone software, algorithm, or computational method was developed for this study. This repository contains custom R analysis scripts used to implement established single-nucleus and spatial transcriptomic workflows, perform statistical analyses, and generate manuscript-associated results.

The scripts are study-specific analysis records rather than an installable R package or a one-command software pipeline. They operate on processed expression matrices, Seurat objects, spatial coordinates, sample metadata, and intermediate files generated during the study.

## Repository contents

| Script | Main purpose |
|---|---|
| `SN_newData_subset—updated.R` | Import, quality control, normalization, clustering, annotation, and cell-type-specific analysis of mouse single-nucleus RNA-seq data. |
| `Aging_score.R` | Identification of age-associated transcriptional changes in white-matter Visium data using cell/spot-level and pseudobulk analyses; generation and summarization of ageing-associated scores. |
| `LPS_merfish.R` | Quality control, clustering, annotation, and differential-expression/module-score analyses of MERFISH data from LPS-treated and control samples. |
| `oligo_niche_analysis.R` | Spatial nearest-neighbour analysis of oligodendrocyte states and their local cellular niches in MERFISH data. |
| `Monkey_4analysis.R` | Processing, quality control, clustering, spatial-domain analysis, and age comparison of non-human-primate spatial transcriptomic data. |

## Data availability

The repository does not contain the full study datasets because the analyses use large, study-specific single-nucleus and spatial transcriptomic objects.

Raw and processed data are available as described in the **Data Availability** section of the manuscript:

- Mouse single-nucleus RNA-seq/MERFISH data/Visium data: `[GSE329620]`

The input paths in the scripts refer to the authors’ original analysis environment and must be replaced with local paths to the corresponding downloaded or generated files.

## System requirements

The analyses were performed in R on macOS compute (Mac Studio M1Max).

- R version: `4.5.0 (2025-04-11)`
- Operating system: `13.2.1 (22D68)`

The principal R package dependencies are:

### CRAN packages

- `Seurat`
- `ggplot2`
- `data.table`
- `pheatmap`
- `reshape2`
- `dplyr`
- `tidyr`
- `RANN`
- `scales`

### Bioconductor packages

- `MAST`
- `DESeq2`
- `SpatialExperiment`
- `Banksy`

all other package informations were involved in sessionInfo.txt

## Expected outputs

Depending on the script, the analyses produce:

- quality-control summaries;
- normalized and annotated Seurat objects (`.rds`);
- cluster and cell-state assignments;
- differential-expression result tables;
- age-associated gene lists and scores;
- spatial-neighbourhood composition tables;
- UMAP, heat-map, spatial-map, and statistical summary figures.

Output file names and save locations are specified within the individual scripts and should be changed to match the user’s local directory structure.

## Demonstration data

A separate demonstration dataset is not included. The repository contains analysis scripts rather than standalone software, and the analyses depend on large, study-specific expression objects, spatial coordinates, and metadata. The scripts are therefore intended to be used with the study data described in the manuscript and its Data Availability statement.

## Reproducibility notes

- No novel algorithm was introduced in this study.
- Analysis parameters, filtering thresholds, clustering resolutions, gene sets, and statistical cut-offs are recorded directly in the scripts.

## Contact

For questions regarding the analysis code, please contact:

**Ruoqing Feng**  
Email: fengruoqing123@gmail.com
