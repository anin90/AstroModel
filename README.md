# AstroModel
This pipeline implements the method described in 
[IDENTIFYING METABOLIC FLUXES ASSOCIATED WITH LITHIUM RESPONSE IN BIPOLAR DISORDER](https://anin90.github.io/).

## Requirements
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

## Input files and directory tree
Directory names can be changed in the [configuration](#configuration) file
````````````
1.matrix2model/	#extract draft models using MEMs

    1.Zhang/  (#v12 implements the method described in manuscript)
		* v12/abs/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_abs.m (MEMs_fpkm_abs) #code
			* matrix2models_abs_v12.mat (Models_fpkm_abs) #output
		* v12/norm_t1/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_norm_t1.m (MEMs_fpkm_abs) #code
			* matrix2models_norm_t1_v12.mat (Models_fpkm_norm_t1) #output	
		* v12/norm_t2/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_norm_t2.m (MEMs_fpkm_abs) #code
			* matrix2models_norm_t2_v12.mat (Models_fpkm_norm_t1) #output
			
    2.Vadodaria/
    
    3.Koskuvi/

2.expandModel/	#expand draft models using literature driven- and nutrient media constraints
    1.Zhang/
    2.Vadodaria/
    3.Koskuvi/

3.analyzeModel/	#analyze models using FVA & MTA
	1.Vadodaria/
	2.Koskuvi/

4.modelComparison/	#compare models with predecessors  

5.generateFigures/	#reproduce figures for manuscript  

````````````

## Configuration

## Running instructions
   * Inputs preparation
   * Run pipeline

## Results

