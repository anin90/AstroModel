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
## Input 'filename': expressionMatrix
## Function: matrix2models_xxx()
## Output: draft models generated using MEMs (n=4).
## (replace 'xxx' with 'abs' (or) 'norm_t1' (or) 'norm_t2', 
##  to generate models for the respective fpkm thresholding).

## Zhang
## Primary_Ctrl
## Input 'filename': 'GSE73721_HMA_CTX.mat'
	cd  1.matrix2model/1.Zhang/v12/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx(filename)
		
## Vadodaria
## iPS_Ctrl
## Input 'filename': 'Vadodaria_Control_Untreated.mat'
	cd  1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx_vadodaria(filename)

## Vadodaria
## iPS_BD
## Input 'filename': 'Vadodaria_BD_Untreated.mat'
	cd  1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx_vadodaria(filename)

## Vadodaria
## iPS_BD_R
## Input 'filename': 'Vadodaria_BD_Untreated_Responder.mat'
	cd  1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx_vadodaria(filename)

## Vadodaria
## iPS_BD_NR
## Input 'filename': 'Vadodaria_BD_Untreated_NonResponder.mat'
	cd  1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/xxx/
	[iMAT_model_TP, GIMME_model_TP, MBA_model_TP, FastCore_model_TP] = matrix2models_xxx_vadodaria(filename)

## Koskuvi
## Healthy Control
## Input 'filename': 'Koskuvi_Control.mat'
	cd  1.matrix2model/3.Koskuvi/1.Control/v1/xxx/
	[iMAT_model_TP] = matrix2models_xxx_koskuvi(filename)

## Koskuvi
## Healthy Twin
## Input 'filename': 'Koskuvi_HT.mat'
	cd  1.matrix2model/3.Koskuvi/2.HT/v1/xxx/
	[iMAT_model_TP] = matrix2models_xxx_koskuvi(filename)

## Koskuvi
## Schiz Twin
## Input 'filename': 'Koskuvi_ST.mat'
	cd  1.matrix2model/3.Koskuvi/3.ST/v1/xxx/
	[iMAT_model_TP] = matrix2models_xxx_koskuvi(filename)

```
   * Generate model statistics
```matlab
## Generate model statistics for all draft models:
## script: generateModelStatsMatrix.m
## Output dir: 1.matrix2model/modelStatsMatFiles/
## Output filename: modelStatsMatSol_YYY.csv (YYY: phenotype)
## Output colname: 'model_ID', 'modelMets', 'modelRxns', 'modelGenes', 'fluxInconsistentRxns',
## Output colname: 'coreRxns', 'overlapCoreRxns', 'overlapCoreRxnsPrct', 'astroModelLewisRxns',
## Output colname: 'overlapLewis', 'overlapLewisPrct'
	
	cd 1.matrix2model/
	run generateModelStatsMatrix.m
```

### 2. expandModel/ - <ins>Expansion of draft models</ins>.
   * Test
   * Test
   * Test

### 3. analyzeModel/ - <ins>Identifying disrupted reactions & subSystems in BD.</ins>

   * Run FVA & MTA to identify reactions disrupted in "BD-lumped", "BD-Responders" and "BD-NonResponders".
```matlab
# Print 'rxnID' & 'subSystems' for each model:
## script: annotateRxnSubsystems.m
## Output dir: 3.analyzeModel/Annotations/AnnotateRxnSubsystems/
## Output filename: Rxns_model_TP_xxx.csv
## where 'model' ~ Primary_Ctrl/iPS_Ctrl/iPS_BD/iPS_BD_R/iPS_BD_NR
## where 'xxx' ~ abs/norm_t1/norm_t2
## Output colname: 'Rxn', 'SubSystem'
	
	cd 3.analyzeModel/Annotations/AnnotateRxnSubsystems/
	run annotateRxnSubsystems.m
	
#FVA:
## Run FVA & identify rxns with 
## FSr >1.5 & <0.8 between iPS-Ctrl & Primary_Ctrl models, and those
## Unchanged between iPS-Ctrl & Primary_Ctrl model,
## script: analyzeCtrl.m
## Output dir: 3.analyzeModel/1.Vadodaria/FSr_Ctrl/
## Output filename: iAstro_FluxDiff_iPSCtrl_vs_Primary.mat
## Output filename: FSR_iAstro_model-1_TP_xxx_vs_iAstro_model-2_xxx.csv
## Output filename: UnChanged_iAstro_model-1_TP_xxx_vs_iAstro_model-2_TP_xxx.csv
## where 'model-1' ~ iPS-Ctrl
## where 'model-2' ~ Primary_Ctrl
## where 'xxx' ~ abs/norm_t1/norm_t2
## Output colname: 'm' (i.e., 'rxnID')
	
	cd 3.analyzeModel/1.Vadodaria/FSr_Ctrl/
	run analyzeCtrl.m

## Run FVA & identify rxns with FSr >1.5 & <0.8 
## between IPS-Ctrl & (IPS-Ctrl-BD; IPS-Ctrl-BD_R; BD_NR):
## script: analyzeBD.m
## Output dir: 3.analyzeModel/1.Vadodaria/FSr_BD/
## Output filename: iAstro_FluxDiff_BD.mat
## Output filename: FSR_iAstro_model-1_TP_xxx_vs_iAstro_model-2_TP_xxx.csv
## where 'model-1' ~ iPS-Ctrl
## where 'model-2' ~ iPS-BD/iPS-BD_R/iPS-BD_NR
## where 'xxx' ~ abs/norm_t1/norm_t2
## Output colname: 'm' (i.e., 'rxnID')
	
	cd 3.analyzeModel/1.Vadodaria/FSr_BD/
	run analyzeBD.m
	
#MTA:
## Run MTA twice by swapping the source and the target states. 
## E.g., in order to identify the reactions relavant to BD, we 
## identified the reactions 'whose knockout' transformed 
## i) “iPS-Control” to “iPS-BD”, and ii) “iPS-BD” to “iPS-Control”. 
## For either runs, the top 20% predictions were first identified, 
## and subsequently their union set were considered 
## for downstream analysis.
## script: runMTA.m
## Output dir: 3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_xxx/
## Output-1 filename: Vadodaria_x_to_y_a.csv (ignore files with '_b')
## Output-2 filename: Vadodaria_y_to_x_a.csv (ignore files with '_b')
## where 'x' ~ Ctrl and 'y' ~ BD/BD_R/BD_NR
## Output colname: 'del_rxnID_KO', 'mta_score', 'alt_score', 'subSystem', 
## Output colname: 'GPR', 'MetabolicUnits', 'Localization'

	cd 3.analyzeModel/1.Vadodaria/MTA_BD/
	run runMTA.m
 ```

   * Filtering reactions relavant to phenotype-of-interest and Reaction-set enrichment analysis (RSEA). 
```matlab

#FVA:
## script: sliceImportantDisruptions_xxx.R
## script: annotateImportantDisruptions_xxx.m
## script: identifyFdrSignificantDisruptions_xxx.R
## script: annotateFinalTable_xxx.m
## Output dir: 3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_xxx/
## Output-1 filename: *_rxns_fdr_xxx.csv
## Output-2 filename: *_subSystem_fdr_xxx.csv
## where * ~ bd_lumped/bd_r/bd_nr
## where 'xxx' ~ abs/norm_t1/norm_t2
## Output-1 colname: 'rxnList', 'subSystem', 'GPR', 'Fluxspan_a', 'Fluxspan_b'
## Output-1 colname: 'FluxSpanRatio', 'Flux', 'MetabolicUnits', 'Localization', 'RxnFormula'
## Output-2 colname: 'subSystem', 'n.x', 'n.y', 'p.val', 'p.val.fdr'
 	
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
