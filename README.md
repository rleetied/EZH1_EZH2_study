# EZH1/2 project: processing and analysis code summaries

This directory contains publication-oriented summaries and the associated analysis code used for final Figures 1-4. The files are intended for archival on GitHub and Zenodo. They preserve the software, parameters, thresholds, analytical sequence, and the relevant implementation code while excluding unrelated exploratory work.

Each assay is organized in two layers:

1. `*_processing.sh` and `*_analysis.R` provide a concise, readable overview of the complete workflow.
2. `processing_code/` and `analysis_code/` contain numbered copies of the original project scripts that implement the corresponding steps and final figure analyses.

The overview scripts should therefore be read first; the numbered code modules provide the underlying commands, statistical models, classifications, and plotting code.

Each assay has separate processing and analysis components:

| Assay | Processing | Analysis | Final figures |
|---|---|---|---|
| RNA-seq | `RNA-seq/RNAseq_processing.sh` | `RNA-seq/RNAseq_analysis.R` | Figure 1G-J |
| ATAC-seq | `ATAC-seq/ATACseq_processing.sh` | `ATAC-seq/ATACseq_analysis.R` | Figures 3F-I and 4A-D |
| ChIP-seq | `ChIP-seq/ChIPseq_processing.sh` | `ChIP-seq/ChIPseq_analysis.R` | Figures 3 and 4E-F |
| ChIP-EM-seq | `EM-seq/EMseq_processing.sh` | `EM-seq/EMseq_analysis.R` | Figure 4F-G |
| Multi-omic integration | `Integration/Integration_processing.sh` | `Integration/Integration_analysis.R` | Figure 4C-J |

Paths in the overview scripts are represented by variables such as `${PROJECT_ROOT}`, `${FASTQ_DIR}`, and `PROJECT_ROOT`. Users should replace these with local paths. Commands containing `<...>` indicate dataset-specific inputs. Some implementation modules retain historical absolute paths and intermediate filenames to document the exact analysis that was performed; these are not expected to execute unchanged on another system.

## Code map

### RNA-seq

- `processing_code/01_trim_galore_generator.sh`: adapter and quality trimming.
- `processing_code/02_STAR_two_pass_generator.sh`: STAR two-pass alignment and gene-level counts.
- `analysis_code/01_differential_expression.R`: primary RKO/SW620 expression processing and contrasts.
- `analysis_code/02_PCA_and_TD_VD_distance.R`: shared-scale PCA and treatment-distance analysis.
- `analysis_code/03_VEH_anchored_response_boxplots.R`: line-specific treatment-response summaries.
- `analysis_code/04_top1000_log2FC_heatmaps.R`: top-variable-gene heatmaps.
- `analysis_code/05_feature_cutoff_sensitivity.R`: 1,000/2,000/5,000/10,000-gene sensitivity analysis.
- `analysis_code/06_selected_pathways.R`: selected pathway comparison.
- `analysis_code/07_EZH1_EZH2_transcripts.R`: EZH1/EZH2 transcript abundance.

### ATAC-seq

- `processing_code/01_ATAC_upstream_processing.sh`: trimming, alignment, duplicate handling, filtering, peak calling, consensus regions, and counting.
- `analysis_code/01_ATAC_QC_PCA.R`: count import, normalization, QC, and PCA.
- `analysis_code/02_differential_accessibility.R`: DESeq2 contrasts and differential-accessibility calls.
- `analysis_code/03_TD_VD_overlap_20pct.R`: TD/VD overlap and the VD >20% classification.
- `analysis_code/04_magnitude_classes_H3K27me1_distance.R`: comparable and VD-favored classes relative to H3K27me1.
- `analysis_code/05_DAC_OR_and_distance.R`: DAC comparisons, odds ratios, and cumulative distance curves.

### ChIP-seq

- `processing_code/01_build_and_alignment.sh`: upstream ChIP-seq processing retained from the production workflow.
- `processing_code/02_computeMatrix_H3K27me1.sh` and `03_computeMatrix_H3K27me3.sh`: deepTools matrices.
- `analysis_code/01_siQ_H3K27me1_analysis.R`: quantitative H3K27me1 processing and visualization.
- `analysis_code/02_H3K27me1_gain_loss_bars.R`: gain/loss region number and genomic coverage.
- `analysis_code/03_ATAC_H3K27me1_distance.R`: accessibility proximity to H3K27me1 dynamics.

### ChIP-EM-seq

- `processing_code/01_trim_galore_generator.sh`: trimming.
- `processing_code/02_BISCUIT_alignment.sh`: BISCUIT alignment and duplicate handling.
- `processing_code/03_BISCUIT_pileup.sh`: methylation pileup generation.
- `processing_code/04_merge_CpG.sh`: strand-aware CpG merging.
- `analysis_code/01_fragment_methylation_thresholds.R`: fragment-level methylation classification.
- `analysis_code/02_waterfall_500bp_50pct_status.R`: position-resolved methylation heatmaps and 50% status bar.

### Multi-omic integration

- `processing_code/01_computeMatrix_ATAC.sh` through `04_computeMatrix_acEMseq.sh`: signal matrices around selected regions.
- `processing_code/05_HOMER_annotation.sh`: genomic annotation and motif-analysis input generation.
- `analysis_code/01_local_rGREAT.R`: local GREAT regulatory-domain associations.
- `analysis_code/02_link_GREAT_to_RNA_priority.R`: GREAT-to-RNA linking and duplicate-gene priority rules.
- `analysis_code/03_barrier_reclassification_scatter.R`: chromatin-barrier classification and integrated scatterplots.
- `analysis_code/04_export_heatmap_annotations.R`: gene/barrier annotation export.
- `analysis_code/05_export_GREAT_expression_table.R`: final integrated gene table.
- `analysis_code/06_ATAC_RNA_scatterplots.R`: accessibility-versus-expression figures.

