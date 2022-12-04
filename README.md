# AstroModel
This pipeline implements the method described in 
[IDENTIFYING METABOLIC FLUXES ASSOCIATED WITH LITHIUM RESPONSE IN BIPOLAR DISORDER](file:///home/anirudh/Downloads/3.-Identifying-Metabolic_Anirudh-S.-Chellappa-Abstract_GBDS_2022.pdf), 
This pipeline runs on a Linux machine only.

## Requirements
1. MATLAB (R2015b+) packages:
   * COBRA Toolbox (v3.0)
   * Gurobi optimizer ()
   * [Metabolic Transformation Algorithm (MTA)](https://github.com/ImNotaGit/MTA)
   * matplotlib (3.3.3+)
2. R (3.6.3+) packages:
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
5. ~50 GB space, and more space will be required depending on the number of models that are built

## Input files and directory tree
Directory names can be changed in the [configuration](#configuration) file
```
1.matrix2model/					#extract draft models using MEMs 
    1.Zhang/
    2.Vadodaria/
    3.Koskuvi/
	modelStatsMatFiles

2.expandModel/					#expand draft models using literature driven- and nutrient media constraints  
    1.Zhang/
    2.Vadodaria/
    3.Koskuvi/

3.analyzeModel/					#analyze models using FVA & MTA  
    1.Zhang/
    2.Vadodaria/
    3.Koskuvi/

4.modelComparison/			#compare models with predecessors  

5.generateFigures/				#reproduce figures for manuscript  

```

## Configuration

## Running instructions

### Inputs preparation

### Run pipeline

## Results

Footer
