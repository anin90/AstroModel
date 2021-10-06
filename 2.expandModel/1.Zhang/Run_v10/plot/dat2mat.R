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

pdf("plot_fsr_intersect.pdf")

###################################################
#rnxs to matrix - import model constraints - v10_TP
###################################################

Rat_Seq_Rxns_Index = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run_v10/plot/FSr_rxns/v10/model_TP/Akkouh_v10_TP_Input.csv");
Rat_Array_Rxns_Index = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run_v10/plot/FSr_rxns/v10/model_TP/GSE66276_Li_v10_TP_Input.csv");
Mouse_Array_Rxns_Index = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run_v10/plot/FSr_rxns/v10/model_TP/GSE132397_Li_v10_TP_Input.csv");

##########################################
# rnxs to matrix: import output - v10_TP
##########################################

Rat_Seq_Rxns_Secondary = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run_v10/plot/FSr_rxns/v10/model_TP/Akkouh_v10_TP.csv");
Rat_Array_Rxns_Secondary = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run_v10/plot/FSr_rxns/v10/model_TP/GSE66276_Li_v10_TP.csv");
Mouse_Array_Rxns_Secondary = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run_v10/plot/FSr_rxns/v10/model_TP/GSE132397_Li_v10_TP.csv");

###################################
#rnxs to matrix - wrangle - v10_TP
###################################

mat = lst(Rat_Seq_Rxns_Secondary, Rat_Array_Rxns_Secondary, Mouse_Array_Rxns_Secondary,
					Rat_Seq_Rxns_Index, Rat_Array_Rxns_Index, Mouse_Array_Rxns_Index) %>% 
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)
  
mat = t(mat)

write.table(mat, "fsr_v10_TP.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

##################################
#rnxs to matrix - UpSetR - v10_TP
##################################

mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run_v10/plot/fsr_v10_TP.csv",
		header = T, sep = "\t")

p = upset(mat,
	order.by = c("degree","freq"), #query.legend = "bottom",
	matrix.color="black", sets = c("Rat_Seq_Rxns_Secondary", "Rat_Array_Rxns_Secondary", "Mouse_Array_Rxns_Secondary",
	"Rat_Seq_Rxns_Index", "Rat_Array_Rxns_Index", "Mouse_Array_Rxns_Index"), keep.order = TRUE,
	sets.bar.color=c("dark grey","dark grey","dark grey","light blue","light blue","light blue"))
p

svg(filename="fsr_v10_TP.svg")
p 
dev.off()




