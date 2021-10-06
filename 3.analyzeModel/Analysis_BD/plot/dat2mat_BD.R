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

pdf("FSR_UpSet_BD.pdf")

##########################################
# rnxs to matrix: import output
##########################################

iPS_Ctrl_TP_vs_iPS_BD_TP = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_TP.csv");
iPS_Ctrl_TP_vs_iPS_BD_R_TP = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_R_TP.csv");
iPS_Ctrl_TP_vs_iPS_BD_NR_TP = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_NR_TP.csv");
FSR_iPS_Ctrl_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_Ctrl_TP_Akkouh_Li.csv");
FSR_iPS_BD_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_TP_Akkouh_Li.csv");
FSR_iPS_BD_R_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_R_TP_Akkouh_Li.csv");
FSR_iPS_BD_NR_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_NR_TP_Akkouh_Li.csv");
RxnInput_iPS_Ctrl_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_Ctrl_TP_Akkouh_Li.csv");
RxnInput_iPS_BD_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_TP_Akkouh_Li.csv");
RxnInput_iPS_BD_R_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_R_TP_Akkouh_Li.csv");
RxnInput_iPS_BD_NR_TP_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_NR_TP_Akkouh_Li.csv");
UnChanged_Primary_TP_vs_iPS_Ctrl_TP = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Ctrl/UnChanged_iAstro_Primary_TP_vs_iAstro_iPS_Ctrl_TP.csv");

iPS_Ctrl_T_vs_iPS_BD_T = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_T_vs_iAstro_iPS_BD_T.csv");
iPS_Ctrl_T_vs_iPS_BD_R_T = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_T_vs_iAstro_iPS_BD_R_T.csv");
iPS_Ctrl_T_vs_iPS_BD_NR_T = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_T_vs_iAstro_iPS_BD_NR_T.csv");
FSR_iPS_Ctrl_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_Ctrl_T_Akkouh_Li.csv");
FSR_iPS_BD_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_T_Akkouh_Li.csv");
FSR_iPS_BD_R_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_R_T_Akkouh_Li.csv");
FSR_iPS_BD_NR_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_fluxdiff/FSR_iAstro_iPS_BD_NR_T_Akkouh_Li.csv");
RxnInput_iPS_Ctrl_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_Ctrl_T_Akkouh_Li.csv");
RxnInput_iPS_BD_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_T_Akkouh_Li.csv");
RxnInput_iPS_BD_R_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_R_T_Akkouh_Li.csv");
RxnInput_iPS_BD_NR_T_Akkouh = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Rxn_input/RxnInput_iAstro_iPS_BD_NR_T_Akkouh_Li.csv");
UnChanged_Primary_T_vs_iPS_Ctrl_T = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Ctrl/UnChanged_iAstro_Primary_T_vs_iAstro_iPS_Ctrl_T.csv");

###################################
# rnxs to matrix - wrangle
###################################

mat = lst(iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP,
			iPS_Ctrl_T_vs_iPS_BD_T, iPS_Ctrl_T_vs_iPS_BD_R_T, iPS_Ctrl_T_vs_iPS_BD_NR_T,
			
			FSR_iPS_Ctrl_TP_Akkouh, FSR_iPS_BD_TP_Akkouh, FSR_iPS_BD_R_TP_Akkouh, FSR_iPS_BD_NR_TP_Akkouh,
			FSR_iPS_Ctrl_T_Akkouh, FSR_iPS_BD_T_Akkouh, FSR_iPS_BD_R_T_Akkouh, FSR_iPS_BD_NR_T_Akkouh,
			
			RxnInput_iPS_Ctrl_TP_Akkouh, RxnInput_iPS_BD_TP_Akkouh, RxnInput_iPS_BD_R_TP_Akkouh, RxnInput_iPS_BD_NR_TP_Akkouh,
			RxnInput_iPS_Ctrl_T_Akkouh, RxnInput_iPS_BD_T_Akkouh, RxnInput_iPS_BD_R_T_Akkouh, RxnInput_iPS_BD_NR_T_Akkouh, 
			
			UnChanged_Primary_TP_vs_iPS_Ctrl_TP, UnChanged_Primary_T_vs_iPS_Ctrl_T) %>% 
			
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)
  
mat = t(mat)

write.table(mat, "FSR_Table_BD.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/FSR_Table_BD.csv",
		header = T, sep = "\t")

##################################
# rnxs to matrix - UpSetR (TP)
##################################

p = upset(mat,
	order.by = c("degree","freq"), #query.legend = "bottom",
	matrix.color="black", sets = c("UnChanged_Primary_TP_vs_iPS_Ctrl_TP", "iPS_Ctrl_TP_vs_iPS_BD_TP", 
									"iPS_Ctrl_TP_vs_iPS_BD_R_TP", "iPS_Ctrl_TP_vs_iPS_BD_NR_TP"), 
	keep.order = TRUE, 

	sets.bar.color=c("dark grey","dark grey","dark grey","dark blue"), nintersects = 60,
	queries = list(list(query = intersects, params = list("iPS_Ctrl_TP_vs_iPS_BD_NR_TP", "UnChanged_Primary_TP_vs_iPS_Ctrl_TP"), color = "red", active = T), 
	list(query = intersects, params = list("iPS_Ctrl_TP_vs_iPS_BD_R_TP", "UnChanged_Primary_TP_vs_iPS_Ctrl_TP"), color = "blue", active = T)))
p

svg(filename="FSR_UpSet_BD_TP.svg")
p 
dev.off()

################
# rnxs slice
################

#bd_17
x = list(UnChanged_Primary_TP_vs_iPS_Ctrl_TP, iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
xInt = Reduce(intersect, x)
y = iPS_Ctrl_TP_vs_iPS_BD_NR_TP
bd_17 = setdiff(xInt,y)
length(bd_17)
write.table(bd_17, "bd_17.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_33
x = list(UnChanged_Primary_TP_vs_iPS_Ctrl_TP, iPS_Ctrl_TP_vs_iPS_BD_TP,
			iPS_Ctrl_TP_vs_iPS_BD_R_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
bd_33 = Reduce(intersect, x)		
length(bd_33)
write.table(bd_33, "bd_33.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_59
x = list(iPS_Ctrl_TP_vs_iPS_BD_NR_TP, iPS_Ctrl_TP_vs_iPS_BD_TP, UnChanged_Primary_TP_vs_iPS_Ctrl_TP)
xInt = Reduce(intersect, x)
y = iPS_Ctrl_TP_vs_iPS_BD_R_TP
bd_59 = setdiff(xInt,y)
length(bd_59)
write.table(bd_59, "bd_59.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_103
x = list(UnChanged_Primary_TP_vs_iPS_Ctrl_TP, iPS_Ctrl_TP_vs_iPS_BD_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_R_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
yInt = Reduce(union,y)
bd_103 = setdiff(xInt,yInt)
length(bd_103)
write.table(bd_103, "bd_103.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_r_92
x = list(UnChanged_Primary_TP_vs_iPS_Ctrl_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
yInt = Reduce(union,y)
bd_r_92 = setdiff(xInt,yInt)
length(bd_r_92)
write.table(bd_r_92, "bd_r_92.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_nr_670
x = list(UnChanged_Primary_TP_vs_iPS_Ctrl_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
yInt = Reduce(union,y)
bd_nr_670 = setdiff(xInt,yInt)
length(bd_nr_670)
write.table(bd_nr_670, "bd_nr_670.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)








