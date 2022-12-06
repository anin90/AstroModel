## Astrocyte metabolic modeling - Identifying metabolic fluxes associated with lithium response in bipolar disorder.
This pipeline implements the method described in [this manuscript](https://anin90.github.io/).

## Requirements:
1. MATLAB (R2015b+):
   * [Gurobi optimizer](https://www.gurobi.com/downloads/licenses/) (obtain license)
	   ```shell
	   cd /opt/gurobi9.5.2_linux64/gurobi952
	   ./grbgetkey "licenseID"
		```
   * [COBRA Toolbox (v3.0)](https://opencobra.github.io/cobratoolbox/stable/installation.html)
   * [Metabolic Transformation Algorithm (MTA)](https://github.com/ImNotaGit/MTA)
2. R (3.6.3+):
   * pacman ()
	   ```r
	   ## check for missing required packages, install them.
	   pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
			readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
			ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
			mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table)
		```						
4. NGS analysis tools: FastQC (), Cutadapt (), HISAT2 (), Samtools (), Cufflinks ()

5. ~5 GB space, and more space will be required depending on the number of models that are built

## Running instructions:
### Clone this repository:
```shell
git clone https://github.com/anin90/AstroModel
cd AstroModel/
```
### 0. prepareExpressionMatrix/ - <ins>Pre-processing and QC of transcriptomes</ins>.
   * Test
   * Test
   * Test
   
### 1. matrix2model/ - <ins>Integration of phenotype-specific transcriptomes with Recon3D</ins>.
   * Extract draft models using "MEMs" (iMAT, GIMME, MBA, FastCore)
```matlab
## (replace 'xxx' with 'abs' (or) 'norm_t1' (or) 'norm_t2', 
##  to generate results for the respective models)	

## Zhang
	## Primary_astrocyte
	## Input 'filename': expressionMatrix (i.e., GSE73721_HMA_CTX.mat)
	cd  1.matrix2model/1.Zhang/v12/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx(filename)
		
## Vadodaria
	## iPS_Ctrl
	## Input 'filename': expressionMatrix (i.e., Vadodaria_Control_Untreated.mat)
	cd  1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx_vadodaria(filename)

	## iPS_BD
	## Input 'filename': expressionMatrix (i.e., Vadodaria_BD_Untreated.mat)
	cd  1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx_vadodaria(filename)

	## iPS_BD_R
	## Input 'filename': expressionMatrix (i.e., Vadodaria_BD_Untreated_Responder.mat)
	cd  1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx_vadodaria(filename)

	## iPS_BD_NR
	## Input 'filename': expressionMatrix (i.e., Vadodaria_BD_Untreated_NonResponder.mat)
	cd  1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx_vadodaria(filename)
	
## Koskuvi
	## Healthy Control
	## Input 'filename': expressionMatrix (i.e., Koskuvi_Control.mat)
	cd  1.matrix2model/3.Koskuvi/1.Control/v1/xxx/
	[iMAT_model_TP] = matrix2models_xxx_koskuvi(filename)

	## Healthy Twin
	## Input 'filename': expressionMatrix (i.e., Koskuvi_HT.mat)
	cd  1.matrix2model/3.Koskuvi/2.HT/v1/xxx/
	[iMAT_model_TP] = matrix2models_xxx_koskuvi(filename)

	## Schiz Twin
	## Input 'filename': expressionMatrix (i.e., Koskuvi_ST.mat)
	cd  1.matrix2model/3.Koskuvi/3.ST/v1/xxx/
	[iMAT_model_TP] = matrix2models_xxx_koskuvi(filename)

```
   * Generate model statistics
```matlab

```

### 2. expandModel/ - <ins>Expansion of draft models</ins>.
   * Test
   * Test
   * Test

### 3. analyzeModel/ - <ins>Identifying disrupted reactions & subSystems in BD.</ins>

   * Run FVA & MTA to identify reactions disrupted in "BD-lumped", "BD-Responders" and "BD-NonResponders".
```matlab
#FVA:
	cd 3.analyzeModel/1.Vadodaria/FSr_Ctrl/
	run analyzeCtrl.m
	
	cd 3.analyzeModel/1.Vadodaria/FSr_BD/
	run analyzeBD.m
	
	cd 3.analyzeModel/Annotations/AnnotateRxnSubsystems/
	run annotateRxnSubsystems.m
#MTA:
	cd 3.analyzeModel/1.Vadodaria/MTA_BD/
	run runMTA.m
 ```

   * Filtering reactions relavant to phenotype-of-interest and Reaction-set enrichment analysis (RSEA). 
```matlab
## (replace 'xxx' with 'abs' (or) 'norm_t1' (or) 'norm_t2', 
##  to generate results for the respective models)	
#FVA:
	cd 3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/
	Rscript sliceImportantDisruptions_xxx.R
	run annotateImportantDisruptions_xxx.m
	Rscript identifyFdrSignificantDisruptions_xxx.R
	run PlotResults/annotateFinalTable_xxx.m
	Rscript plotFinalTable_xxx.R
#MTA:
	cd 3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/
	Rscript analyzeMTAscores_xxx.R
```

   * Identifying disruptions that are significant across modules.
```matlab
	cd 3.analyzeModel/3.DataOverlays/
	run generateDataOverlays.m
	Rscript plotDataOverlays.R
	run annotateDataOverlays.m
	Rscript filterDataOverlays.R 
```

### 4. modelComparison/ - <ins>Compare model statistics</ins>.
   * Test
   * Test
   * Test
   
### 5. generateFigures/ - <ins>Generate figures for manuscript</ins>.
   * Test
   * Test
   * Test
