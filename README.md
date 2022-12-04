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

## Files (data, code and output) and directory tree
````````````
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
0.prepareExpressionMatrix/	#transcriptome processing & QC
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
1.matrix2model/	#extract draft models using MEMs
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    1.Zhang/

		(#phenotype- primary astrocytes)
		(#v12 implements the method described in manuscript)
		
		* v12/abs/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_abs.m (MEMs_fpkm_abs) #code
			* matrix2models_abs_v12.mat (Models_fpkm_abs) #output
		* v12/norm_t1/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_norm_t1.m (MEMs_fpkm_norm_t1) #code
			* matrix2models_norm_t1_v12.mat (Models_fpkm_norm_t1) #output	
		* v12/norm_t2/
			* GSE73721_HMA_CTX.mat (ExpressionMatrix) #data
			* matrix2models_norm_t2.m (MEMs_fpkm_norm_t2) #code
			* matrix2models_norm_t2_v12.mat (Models_fpkm_norm_t1) #output
			
    2.Vadodaria/

		(#phenotype- iPS-derived astrocytes from BD patients and controls)
		
		(#Below mentioned only for "1.Control_Untreated". In similar logic, 
		data & codes for "2.BD_Untreated", "3.BD_Responder_Untreated" 
		and "4.BD_NonResponder_Untreated" are available under "2.Vadodaria/")
				
		(#v3 implements the method described in manuscript)
		
		* 1.Control_Untreated/v3/abs/
			* Vadodaria_Control_Untreated.mat (ExpressionMatrix) #data
			* matrix2models_abs_vadodaria.m (MEMs_fpkm_abs) #code
			* matrix2models_abs_v3.mat (Models_fpkm_abs) #output
		* 1.Control_Untreated/v3/norm_t1/
			* Vadodaria_Control_Untreated.mat (ExpressionMatrix) #data
			* matrix2models_norm_t1_vadodaria.m (MEMs_fpkm_norm_t1) #code
			* matrix2models_norm_t1_v3.mat (Models_fpkm_norm_t1) #output
		* 1.Control_Untreated/v3/norm_t2/
			* Vadodaria_Control_Untreated.mat (ExpressionMatrix) #data
			* matrix2models_norm_t2_vadodaria.m (MEMs_fpkm_norm_t2) #code
			* matrix2models_norm_t2_v3.mat (Models_fpkm_norm_t2) #output
    
    3.Koskuvi/

		(#phenotype- iPS-derived astrocytes from monozygotic twin pairs discordant 
		for schizophrenia and healthy subjects)
		
		(#Below mentioned only for "1.Control"  (healthy controls). In similar logic, 
		data & codes for "2.HT" (healthy-twin) and "3.ST" (schizophrenia-twin) are 
		available under "3.Koskuvi/")

		(#v1 implements the method described in manuscript)

		* 1.Control/v1/abs/
			* Koskuvi_Control.mat (ExpressionMatrix) #data
			* matrix2models_abs_koskuvi.m (MEMs_fpkm_abs) #code
			* matrix2models_abs_v1.mat (Models_fpkm_abs) #output
		* 1.Control/v1/norm_t1/
			* Koskuvi_Control.mat (ExpressionMatrix) #data
			* matrix2models_norm_t1_koskuvi.m (MEMs_fpkm_norm_t1) #code
			* matrix2models_norm_t1_v1.mat (Models_fpkm_norm_t1) #output
		* 1.Control/v1/norm_t2/
			* Koskuvi_Control.mat (ExpressionMatrix) #data
			* matrix2models_norm_t2_koskuvi.m (MEMs_fpkm_norm_t2) #code
			* matrix2models_norm_t2_v1.mat (Models_fpkm_norm_t2) #output

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
2.expandModel/	#expand draft models using literature & nutrient media constraints
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    1.Zhang/
    2.Vadodaria/
    3.Koskuvi/

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
3.analyzeModel/	#analyze models using FVA & MTA
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
4.modelComparison/	#compare models with predecessors  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
5.generateFigures/	#generate/reproduce figures for manuscript
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


````````````

## Configuration

## Running instructions
   * Inputs preparation
   * Run pipeline

## Results

