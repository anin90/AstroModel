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
   * Test
   * Test
   * Test

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
```
## ('xxx' - abs/norm_t1/norm_t2)	
#FVA:
	Rscript 3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/sliceImportantDisruptions_xxx.R
	run 3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/annotateImportantDisruptions_xxx.m
	Rscript 3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/identifyFdrSignificantDisruptions_xxx.R
	run 3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/annotateFinalTable_xxx.m
	Rscript 3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/plotFinalTable_xxx.R
#MTA:
	Rscript 3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/analyzeMTAscores_xxx.R
```

   * Identifying disruptions that are significant across modules.
```
	run 3.analyzeModel/3.DataOverlays/generateDataOverlays.m
	Rscript 3.analyzeModel/3.DataOverlays/plotDataOverlays.R
	run 3.analyzeModel/3.DataOverlays/annotateDataOverlays.m
	Rscript 3.analyzeModel/3.DataOverlays/filterDataOverlays.R 
```

### 4. modelComparison/ - <ins>Compare model statistics</ins>.
   * Test
   * Test
   * Test
   
### 5. generateFigures/ - <ins>Generate figures for manuscript</ins>.
   * Test
   * Test
   * Test
