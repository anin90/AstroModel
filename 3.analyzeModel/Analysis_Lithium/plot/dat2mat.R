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

pdf("FSR_UpSet.pdf")

###################################################
#rnxs to matrix - import model constraints
###################################################

RxnInput_iPS_BD_NR_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_NR_T_Akkouh_Li.csv");
RxnInput_iPS_BD_NR_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_NR_T_GSE132397_Li.csv");
RxnInput_iPS_BD_NR_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_NR_T_GSE66276_Li.csv");
RxnInput_iPS_BD_NR_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_NR_TP_Akkouh_Li.csv");
RxnInput_iPS_BD_NR_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_NR_TP_GSE132397_Li.csv");
RxnInput_iPS_BD_NR_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_NR_TP_GSE66276_Li.csv");
RxnInput_iPS_BD_R_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_R_T_Akkouh_Li.csv");
RxnInput_iPS_BD_R_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_R_T_GSE132397_Li.csv");
RxnInput_iPS_BD_R_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_R_T_GSE66276_Li.csv");
RxnInput_iPS_BD_R_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_R_TP_Akkouh_Li.csv");
RxnInput_iPS_BD_R_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_R_TP_GSE132397_Li.csv");
RxnInput_iPS_BD_R_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_R_TP_GSE66276_Li.csv");
RxnInput_iPS_BD_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_T_Akkouh_Li.csv");
RxnInput_iPS_BD_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_T_GSE132397_Li.csv");
RxnInput_iPS_BD_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_T_GSE66276_Li.csv");
RxnInput_iPS_BD_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_TP_Akkouh_Li.csv");
RxnInput_iPS_BD_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_TP_GSE132397_Li.csv");
RxnInput_iPS_BD_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_TP_GSE66276_Li.csv");
RxnInput_iPS_Ctrl_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_Ctrl_T_Akkouh_Li.csv");
RxnInput_iPS_Ctrl_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_Ctrl_T_GSE132397_Li.csv");
RxnInput_iPS_Ctrl_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_Ctrl_T_GSE66276_Li.csv");
RxnInput_iPS_Ctrl_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_Ctrl_TP_Akkouh_Li.csv");
RxnInput_iPS_Ctrl_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_Ctrl_TP_GSE132397_Li.csv");
RxnInput_iPS_Ctrl_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_Ctrl_TP_GSE66276_Li.csv");
RxnInput_Primary_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_Primary_T_Akkouh_Li.csv");
RxnInput_Primary_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_Primary_T_GSE132397_Li.csv");
RxnInput_Primary_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_Primary_T_GSE66276_Li.csv");
RxnInput_Primary_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_Primary_TP_Akkouh_Li.csv");
RxnInput_Primary_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_Primary_TP_GSE132397_Li.csv");
RxnInput_Primary_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_Primary_TP_GSE66276_Li.csv");

##########################################
# rnxs to matrix: import output
##########################################

FSR_iPS_BD_NR_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_NR_T_Akkouh_Li.csv");
FSR_iPS_BD_NR_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_NR_T_GSE132397_Li.csv");
FSR_iPS_BD_NR_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_NR_T_GSE66276_Li.csv");
FSR_iPS_BD_NR_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_NR_TP_Akkouh_Li.csv");
FSR_iPS_BD_NR_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_NR_TP_GSE132397_Li.csv");
FSR_iPS_BD_NR_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_NR_TP_GSE66276_Li.csv");
FSR_iPS_BD_R_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_R_T_Akkouh_Li.csv");
FSR_iPS_BD_R_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_R_T_GSE132397_Li.csv");
FSR_iPS_BD_R_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_R_T_GSE66276_Li.csv");
FSR_iPS_BD_R_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_R_TP_Akkouh_Li.csv");
FSR_iPS_BD_R_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_R_TP_GSE132397_Li.csv");
FSR_iPS_BD_R_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_R_TP_GSE66276_Li.csv");
FSR_iPS_BD_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_T_Akkouh_Li.csv");
FSR_iPS_BD_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_T_GSE132397_Li.csv");
FSR_iPS_BD_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_T_GSE66276_Li.csv");
FSR_iPS_BD_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_TP_Akkouh_Li.csv");
FSR_iPS_BD_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_TP_GSE132397_Li.csv");
FSR_iPS_BD_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_TP_GSE66276_Li.csv");
FSR_iPS_Ctrl_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_Ctrl_T_Akkouh_Li.csv");
FSR_iPS_Ctrl_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_Ctrl_T_GSE132397_Li.csv");
FSR_iPS_Ctrl_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_Ctrl_T_GSE66276_Li.csv");
FSR_iPS_Ctrl_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_Ctrl_TP_Akkouh_Li.csv");
FSR_iPS_Ctrl_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_Ctrl_TP_GSE132397_Li.csv");
FSR_iPS_Ctrl_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_Ctrl_TP_GSE66276_Li.csv");
FSR_Primary_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_Primary_T_Akkouh_Li.csv");
FSR_Primary_T_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_Primary_T_GSE132397_Li.csv");
FSR_Primary_T_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_Primary_T_GSE66276_Li.csv");
FSR_Primary_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_Primary_TP_Akkouh_Li.csv");
FSR_Primary_TP_GSE132397 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_Primary_TP_GSE132397_Li.csv");
FSR_Primary_TP_GSE66276 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_Primary_TP_GSE66276_Li.csv");

###################################
#rnxs to matrix - wrangle
###################################

mat = lst(RxnInput_Primary_TP_Akkouh, RxnInput_Primary_TP_GSE132397, RxnInput_Primary_TP_GSE66276,
			RxnInput_iPS_Ctrl_TP_Akkouh, RxnInput_iPS_Ctrl_TP_GSE132397, RxnInput_iPS_Ctrl_TP_GSE66276,
			RxnInput_iPS_BD_TP_Akkouh, RxnInput_iPS_BD_TP_GSE132397, RxnInput_iPS_BD_TP_GSE66276,
			RxnInput_iPS_BD_R_TP_Akkouh, RxnInput_iPS_BD_R_TP_GSE132397, RxnInput_iPS_BD_R_TP_GSE66276, 
			RxnInput_iPS_BD_NR_TP_Akkouh, RxnInput_iPS_BD_NR_TP_GSE132397, RxnInput_iPS_BD_NR_TP_GSE66276, 

			FSR_Primary_TP_Akkouh, FSR_Primary_TP_GSE132397, FSR_Primary_TP_GSE66276,
			FSR_iPS_Ctrl_TP_Akkouh, FSR_iPS_Ctrl_TP_GSE132397, FSR_iPS_Ctrl_TP_GSE66276,
			FSR_iPS_BD_TP_Akkouh, FSR_iPS_BD_TP_GSE132397, FSR_iPS_BD_TP_GSE66276,
			FSR_iPS_BD_R_TP_Akkouh, FSR_iPS_BD_R_TP_GSE132397, FSR_iPS_BD_R_TP_GSE66276, 
			FSR_iPS_BD_NR_TP_Akkouh, FSR_iPS_BD_NR_TP_GSE132397, FSR_iPS_BD_NR_TP_GSE66276, 
			
			
			RxnInput_Primary_T_Akkouh, RxnInput_Primary_T_GSE132397, RxnInput_Primary_T_GSE66276,
			RxnInput_iPS_Ctrl_T_Akkouh, RxnInput_iPS_Ctrl_T_GSE132397, RxnInput_iPS_Ctrl_T_GSE66276,
			RxnInput_iPS_BD_T_Akkouh, RxnInput_iPS_BD_T_GSE132397, RxnInput_iPS_BD_T_GSE66276,
			RxnInput_iPS_BD_R_T_Akkouh, RxnInput_iPS_BD_R_T_GSE132397, RxnInput_iPS_BD_R_T_GSE66276, 
			RxnInput_iPS_BD_NR_T_Akkouh, RxnInput_iPS_BD_NR_T_GSE132397, RxnInput_iPS_BD_NR_T_GSE66276,
			
			FSR_Primary_T_Akkouh, FSR_Primary_T_GSE132397, FSR_Primary_T_GSE66276,
			FSR_iPS_Ctrl_T_Akkouh, FSR_iPS_Ctrl_T_GSE132397, FSR_iPS_Ctrl_T_GSE66276,
			FSR_iPS_BD_T_Akkouh, FSR_iPS_BD_T_GSE132397, FSR_iPS_BD_T_GSE66276,
			FSR_iPS_BD_R_T_Akkouh, FSR_iPS_BD_R_T_GSE132397, FSR_iPS_BD_R_T_GSE66276, 
			FSR_iPS_BD_NR_T_Akkouh, FSR_iPS_BD_NR_T_GSE132397, FSR_iPS_BD_NR_T_GSE66276) %>% 
			
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)
  
mat = t(mat)

write.table(mat, "FSR_Table.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

##################################
#rnxs to matrix - UpSetR (all)
##################################

mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/FSR_Table.csv",
		header = T, sep = "\t")

p = upset(mat,
	order.by = c("degree","freq"), #query.legend = "bottom",
	matrix.color="black", sets = c("FSR_iPS_Ctrl_T_Akkouh", "FSR_iPS_Ctrl_T_GSE132397", "FSR_iPS_Ctrl_T_GSE66276",
	"FSR_iPS_BD_T_Akkouh", "FSR_iPS_BD_T_GSE132397", "FSR_iPS_BD_T_GSE66276",
	"FSR_iPS_BD_R_T_Akkouh", "FSR_iPS_BD_R_T_GSE132397", "FSR_iPS_BD_R_T_GSE66276", 
	"FSR_iPS_BD_NR_T_Akkouh", "FSR_iPS_BD_NR_T_GSE132397", "FSR_iPS_BD_NR_T_GSE66276"), keep.order = TRUE,

	sets.bar.color=c("dark grey","dark grey","dark grey",
	"dark grey","dark grey","dark grey",
	"dark grey","dark grey","dark grey",
	"dark grey","dark grey","dark grey"), point.size = 0.9, line.size = 0.9, nintersects = 60)
p

svg(filename="FSR_UpSet.svg")
p 
dev.off()

##################################
#rnxs to matrix - UpSetR (rm_GSE66276)
##################################

mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/FSR_Table.csv",
		header = T, sep = "\t")

p = upset(mat,
	order.by = c("degree","freq"), #query.legend = "bottom",
	matrix.color="black", sets = c("FSR_iPS_Ctrl_T_Akkouh", "FSR_iPS_Ctrl_T_GSE66276",
	"FSR_iPS_BD_T_Akkouh", "FSR_iPS_BD_T_GSE66276",
	"FSR_iPS_BD_R_T_Akkouh", "FSR_iPS_BD_R_T_GSE66276", 
	"FSR_iPS_BD_NR_T_Akkouh", "FSR_iPS_BD_NR_T_GSE66276"), keep.order = TRUE,

	sets.bar.color=c("dark grey","dark grey",
	"dark grey","dark grey",
	"dark grey","dark grey",
	"dark grey","dark grey"), point.size = 0.9, line.size = 0.9, nintersects = 60)
p

svg(filename="FSR_UpSet_rm_GSE66276.svg")
p 
dev.off()

##################################
#rnxs to matrix - UpSetR (empty.intersections = "on")
##################################

mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/FSR_Table.csv",
		header = T, sep = "\t")

p = upset(mat,
	empty.intersections = "on",
	order.by = c("degree","freq"),
	matrix.color="black")
p

svg(filename="FSR_UpSet_keep_empty_intersections.svg")
p 
dev.off()

####
#rnxs intersect
####

#TP
All_TP = list(
			FSR_iPS_Ctrl_TP_Akkouh, FSR_iPS_Ctrl_TP_GSE66276,
			FSR_iPS_BD_TP_Akkouh, FSR_iPS_BD_TP_GSE66276,
			FSR_iPS_BD_R_TP_Akkouh, FSR_iPS_BD_R_TP_GSE66276)
			
NR_TP = list(
							FSR_iPS_BD_NR_TP_Akkouh, FSR_iPS_BD_NR_TP_GSE66276)			
			
All_TP_Intersect = Reduce(intersect, All_TP)					
NR_TP_Intersect = Reduce(intersect, NR_TP)								
			
length(setdiff(All_TP_Intersect, NR_TP_Intersect))	
setdiff(All_TP_Intersect, NR_TP_Intersect)					
			
			
#T
All_T = list(
			FSR_iPS_Ctrl_T_Akkouh, FSR_iPS_Ctrl_T_GSE66276,
			FSR_iPS_BD_T_Akkouh, FSR_iPS_BD_T_GSE66276,
			FSR_iPS_BD_R_T_Akkouh, FSR_iPS_BD_R_T_GSE66276)
			
NR_T = list(
							FSR_iPS_BD_NR_T_Akkouh, FSR_iPS_BD_NR_T_GSE66276)			
			
All_T_Intersect = Reduce(intersect, All_T)					
NR_T_Intersect = Reduce(intersect, NR_T)								
			
length(setdiff(All_T_Intersect, NR_T_Intersect))	
setdiff(All_T_Intersect, NR_T_Intersect)	

