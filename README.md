# RNA-seq DESeq2 Analysis of GSE179354

RNA-seq differential gene expression analysis of the public GEO dataset **GSE179354**
to study cytokine (TNFα, IFNγ, TNFα+IFNγ) responses and the modulatory effects of
fluticasone propionate (FP) in pediatric airway smooth muscle cells.

## Dataset
- GEO accession: **GSE179354**
- Organism: *Homo sapiens*
- Assay: RNA-seq  
Gene-level count data were downloaded from NCBI GEO and processed locally.

## Methods
Lowly expressed genes were filtered prior to analysis.
Differential expression analysis was performed using DESeq2 in R using a main-effects model (~ FP + cytokine) and an interaction model (~ cytokine * FP) to identify cytokine-regulated genes and FP-modulated responses.
Statistical significance was assessed using FDR < 0.05.

## Outputs
- Differential expression results for all contrasts (`.csv`)
- Interaction effect results for FP
- Summary bar plot of up and downregulated genes

## Reproducibility
R session information (R version and package versions) is provided in`Results/sessionInfo.docx`.
