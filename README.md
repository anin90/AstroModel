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

## Files (data, code and output) and directory tree:
````````````
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
0.prepareExpressionMatrix/ #transcriptome processing & QC
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
1.matrix2model/ #extract draft models using "MEMs" (iMAT, GIMME, MBA, FastCore)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   * Zhang/			
   * Vadodaria/
   * Koskuvi/
   * modelStatsMatFiles/

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
2.expandModel/ #expand draft models (only 'iMAT')
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   * Zhang/
   * Vadodaria/
   * Koskuvi/

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
3.analyzeModel/ #analyze models using FVA & MTA
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   * Vadodaria/
   * Koskuvi/
   * DataOverlays/
   * ModelFnTests/
   * Annotations/

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
4.modelComparison/ #compare model statistics  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
5.generateFigures/ #generate figures for manuscript
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


````````````

## Configuration:

## Running instructions:
   * Inputs preparation:
   * Run pipeline:

## Results:

