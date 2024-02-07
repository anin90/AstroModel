library("pheatmap")
library("RColorBrewer")
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
library(splitstackshape)
library(gplots)
library(ggplot2)
library(UpSetR)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(plot.matrix)

pdf("st_tbl_abs/FSR_UpSet_ST_abs.pdf")

##########################################
# rnxs to matrix: import output
##########################################

	iPS_Ctrl_TP_vs_iPS_ST_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/FSR_iAstro_iPS_SCZ_Ctrl_TP_abs_vs_iAstro_iPS_SCZ_ST_TP_abs.csv",header = T, sep = "\t");
	UnChanged_iPS_Ctrl_TP_vs_iPS_HT_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_HT/UnChanged_iAstro_iPS_SCZ_Ctrl_TP_abs_vs_iAstro_iPS_SCZ_HT_TP_abs.csv",header = T, sep = "\t");	
	UnChanged_iPS_Ctrl_TP_vs_Primary_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_Ctrl/UnChanged_iAstro_iPS_SCZ_Ctrl_TP_abs_vs_iAstro_Primary_TP_abs.csv",header = T, sep = "\t");

	Ctrl_vs_ST = iPS_Ctrl_TP_vs_iPS_ST_TP$m
	Ctrl_vs_HT_Unchanged = UnChanged_iPS_Ctrl_TP_vs_iPS_HT_TP$n	
	Ctrl_vs_Primary_Unchanged = UnChanged_iPS_Ctrl_TP_vs_Primary_TP$n

###################################
# rnxs to matrix - wrangle
###################################

	mat = lst(Ctrl_vs_ST, Ctrl_vs_HT_Unchanged, Ctrl_vs_Primary_Unchanged) %>% 
				
	  enframe %>% 
	  unnest %>% 
	  count(name, value) %>% 
	  spread(value, n, fill = 0)
	  
	mat = t(mat)

	write.table(mat, "st_tbl_abs/FSR_Table_ST_abs.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

	mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_abs/FSR_Table_ST_abs.csv",
			header = T, sep = "\t")

##################################
# rnxs to matrix - UpSetR (TP)
##################################

	p = upset(mat,
		order.by = c("degree","freq"), #query.legend = "bottom",
		matrix.color="black", sets = c("Ctrl_vs_Primary_Unchanged", "Ctrl_vs_HT_Unchanged", "Ctrl_vs_ST"), 
		keep.order = TRUE, 
		queries = list(list(query = intersects,  params = list("Ctrl_vs_Primary_Unchanged", "Ctrl_vs_HT_Unchanged", "Ctrl_vs_ST"), 
		color = "red", active = F)), sets.bar.color=c("dark grey", "dark grey", "dark blue"), nintersects = 60)

	p

	svg(filename="st_tbl_abs/FSR_UpSet_ST_TP_abs.svg")
	p 
	dev.off()

################
# rnxs slice
################

	#st_rxns
	x = list(Ctrl_vs_Primary_Unchanged, Ctrl_vs_HT_Unchanged, Ctrl_vs_ST)
	st_rxns = Reduce(intersect, x)		
	length(st_rxns)
	write.table(st_rxns, "st_tbl_abs/st_rxns.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)
