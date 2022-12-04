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
0. prepareExpressionMatrix/ - <ins>Transcriptome processing & QC</ins>:
   * Test
   * Test
   * Test
1. matrix2model/ - <ins>Extract draft models using "MEMs" (iMAT, GIMME, MBA, FastCore)</ins>:
   * Test
   * Test
   * Test
2. expandModel/ - <ins>Expand draft (only iMAT) by adding literature curated and nutrient media constraints</ins>.
   * Test
   * Test
   * Test
3. analyzeModel/ - <ins>Analyze models using FVA & MTA to identify reactions disrupted in BD, BD-Responder & BD-NonResponder</ins>.
   * Test
   * Test
   * Test
4. modelComparison/ - <ins>Compare model statistics</ins>:
   * Test
   * Test
   * Test
5. generateFigures/ - <ins>Generate figures for manuscript</ins>:
   * Test
   * Test
   * Test
