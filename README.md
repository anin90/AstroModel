## Astrocyte metabolic modeling - Identifying metabolic fluxes associated with lithium response in bipolar disorder.
This pipeline implements the method described in [this manuscript](https://anin90.github.io/).

## Requirements:
1. MATLAB (R2015b+):
   * Gurobi optimizer ()
   * COBRA Toolbox (v3.0)
   * [Metabolic Transformation Algorithm (MTA)](https://github.com/ImNotaGit/MTA)
2. R (3.6.3+):
   * pacman ()
   * pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
							readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
							ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
							mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table)
4. NGS analysis tools: 
   * FastQC ()
   * Cutadapt ()
   * HISAT2 ()
   * Samtools ()
   * Cufflinks ()
5. ~5 GB space, and more space will be required depending on the number of models that are built

## Running instructions:
0. prepareExpressionMatrix/ - Transcriptome processing & QC:
   * Test
   * Test
   * Test
1. matrix2model/ #extract draft models using "MEMs" (iMAT, GIMME, MBA, FastCore):
   * Test
   * Test
   * Test
2. expandModel/ #expand draft (only iMAT) by adding literature curated and nutrient media constraints.
   * Test
   * Test
   * Test
3. analyzeModel/ #analyze models using FVA & MTA to identify reactions disrupted in BD, BD_R & BD_NR.
   * Test
   * Test
   * Test
4. modelComparison/ #compare model statistics:
   * Test
   * Test
   * Test
5. generateFigures/ #generate figures for manuscript
   * Test
   * Test
   * Test
