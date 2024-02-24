## Chromatin accessibility atlas of rheumatoid arthritis tissue

Code corresponding to Weinand*, Sakaue* et al., The Chromatin Landscape of Pathogenic Transcriptional Cell States in Rheumatoid Arthritis, Submitted, 2023

Contact: Kathryn Weinand kweinand@fas.harvard.edu 


### Software information:

We used R v3.6.1 for most analyses with the following packages: argparse v2.0.3, aricode v1.0.0, BiocGenerics v0.30.0, class v7.3-17, data.table v1.12.8, dplyr v1.0.2, GenomeInfoDb v1.20.0, GenomicRanges v1.36.1, ggbeeswarm v0.6.0, ggplot2 v3.3.0, ggpubr v0.4.0, ggrastr v0.2.3, ggrepel v0.8.2, ggthemes v4.2.0, gplots v3.0.1.1, gridExtra v2.3, gtools v3.8.2, harmony v1.0, IRanges v2.18.3, irlba v2.3.3, lattice v0.20-41, lme4 v1.1-21, magrittr v1.5, MASS v7.3-51.6, Matrix v1.2-18, Matrix.utils v0.9.7, matrixStats v0.56.0, patchwork v1.1.0.9000, pheatmap v1.0.12, plyr v1.8.6, presto v1.0.0, RANN v2.6.1, RColorBrewer v1.1-2, rcompanion v2.4.1, Rcpp v1.0.4.6, RcppCNPy v0.2.10, repr v1.0.1, reticulate v1.13, Rmisc v1.5.1, ROCR v1.0-7, rstatix v0.7.0, S4Vectors v0.22.1, scales v1.1.1, Seurat v3.2.0, Signac v1.1.0, stringr v1.4.0, symphony v1.0, tibble v3.0.1, tidyr v1.0.3, umap v0.2.3.1, uwot v0.1.8, viridis v0.5.1, viridisLite v0.3.0.

For ArchR analyses, we used R v4.2.0 with the following packages: ArchR v1.0.2, argparse v2.1.6, Biobase v2.56.0, BiocGenerics v0.42.0, Biostrings v2.64.1, BSgenome v1.64.0, BSgenome.Hsapiens.UCSC.hg38 v1.4.4, chromVARmotifs v0.2.0, data.table v1.14.4, GenomeInfoDb v1.32.4, GenomicRanges v1.48.0, ggplot2 v3.3.6, gridExtra v2.3, gtable v0.3.1, gtools v3.9.3, IRanges v2.30.1, JASPAR2016 v1.24.0, JASPAR2018 v1.1.1, JASPAR2020 v0.99.10, magrittr v2.0.3, Matrix v1.5-1, MatrixGenerics v1.8.1, matrixStats v0.62.0, plyr v1.8.7, Rcpp v1.0.9, rhdf5 v2.40.0, rtracklayer v1.56.1, S4Vectors v0.34.0, stringr v1.4.1, SummarizedExperiment v1.26.1, TFBSTools v1.34.0, tidyr v1.2.1, XVector v0.36.0.

We also used python v3.7.3, scrublet v0.2.3, samtools v1.9, bedtools v2.28.0, bedops v2.4.36, GNU Awk 3.1.7, jupyter v4.4.0.

We recommend downloading software via conda.


### Directory information

paper_figures: jupyter notebooks used to create the figures/tables in the submitted manuscript. There are also two R files with the functions used in the jupyter notebooks.

scripts: The scripts used to process the data split up by analysis. Usage information by script can be found by using Rscript {script} -h for R scripts or python {script} -h for python scripts.

