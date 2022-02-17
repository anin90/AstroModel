# http://www.sthda.com/english/wiki/dot-charts-r-base-graphs
# https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
library("pheatmap")
library("RColorBrewer")
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
library(splitstackshape)
library(gplots)
library(ggplot2)
library(ggfortify)
require(reshape2)
library(factoextra)
library(plot.matrix)
library(VennDiagram)
library(ggvenn)
library(plotrix)
library(pheatmap)
library(magrittr)
library(venn)

pdf("bd_MU_tbl.pdf")

###########
# Load data
###########

bd_212_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_212_tbl.csv",
		header = T, sep = "\t")
bd_r_92_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_r_92_tbl.csv",
		header = T, sep = "\t")		
bd_nr_670_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_nr_670_tbl.csv",
		header = T, sep = "\t")

###########
# Plot data
###########
		
# Fig.1: venn rxns
	x <- list(
		BD_Lumped = bd_212_tbl$rxnList, 
		BD_Responder = bd_r_92_tbl$rxnList, 
		BD_NonResponder = bd_nr_670_tbl$rxnList)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlap in rxns disrupted between models")

# Fig.2: venn subsystems
	x <- list(
		BD_Lumped = bd_212_tbl$Var8, 
		BD_Responder = bd_r_92_tbl$Var8, 
		BD_NonResponder = bd_nr_670_tbl$Var8)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlap in MetabolicUnits disrupted between models")				
