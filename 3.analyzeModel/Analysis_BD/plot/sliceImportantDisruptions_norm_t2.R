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

pdf("bd_tbl_norm_t2/FSR_UpSet_BD_norm_t2.pdf")

##########################################
# rnxs to matrix: import output
##########################################

iPS_Ctrl_TP_vs_iPS_BD_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_TP_norm_t2.csv",header = T, sep = "\t");
iPS_Ctrl_TP_vs_iPS_BD_R_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_R_TP_norm_t2.csv",header = T, sep = "\t");
iPS_Ctrl_TP_vs_iPS_BD_NR_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_NR_TP_norm_t2.csv",header = T, sep = "\t");
UnChanged_iPS_Ctrl_TP_vs_Primary_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Ctrl/UnChanged_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_Primary_TP_norm_t2.csv",header = T, sep = "\t");

iPS_Ctrl_TP_vs_iPS_BD_TP = iPS_Ctrl_TP_vs_iPS_BD_TP$m
iPS_Ctrl_TP_vs_iPS_BD_R_TP = iPS_Ctrl_TP_vs_iPS_BD_R_TP$m
iPS_Ctrl_TP_vs_iPS_BD_NR_TP = iPS_Ctrl_TP_vs_iPS_BD_NR_TP$m
UnChanged_iPS_Ctrl_TP_vs_Primary_TP = UnChanged_iPS_Ctrl_TP_vs_Primary_TP$n

###################################
# rnxs to matrix - wrangle
###################################

mat = lst(iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP, 			
			UnChanged_iPS_Ctrl_TP_vs_Primary_TP) %>% 
			
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)
  
mat = t(mat)

write.table(mat, "bd_tbl_norm_t2/FSR_Table_BD_norm_t2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_norm_t2/FSR_Table_BD_norm_t2.csv",
		header = T, sep = "\t")

##################################
# rnxs to matrix - UpSetR (TP)
##################################

p = upset(mat,
	order.by = c("degree","freq"), #query.legend = "bottom",
	matrix.color="black", sets = c("UnChanged_iPS_Ctrl_TP_vs_Primary_TP", "iPS_Ctrl_TP_vs_iPS_BD_TP", 
									"iPS_Ctrl_TP_vs_iPS_BD_R_TP", "iPS_Ctrl_TP_vs_iPS_BD_NR_TP"), 
	keep.order = TRUE, 

	sets.bar.color=c("dark grey","dark grey","dark grey","dark blue"), nintersects = 60,
	queries = list(list(query = intersects, params = list("iPS_Ctrl_TP_vs_iPS_BD_NR_TP", "UnChanged_iPS_Ctrl_TP_vs_Primary_TP"), color = "red", active = T), 
	list(query = intersects, params = list("iPS_Ctrl_TP_vs_iPS_BD_R_TP", "UnChanged_iPS_Ctrl_TP_vs_Primary_TP"), color = "blue", active = T)))
p

svg(filename="bd_tbl_norm_t2/FSR_UpSet_BD_TP_norm_t2.svg")
p 
dev.off()

################
# rnxs slice
################

#bd_79
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
xInt = Reduce(intersect, x)
y = iPS_Ctrl_TP_vs_iPS_BD_NR_TP
bd_79 = setdiff(xInt,y)
length(bd_79)
write.table(bd_79, "bd_tbl_norm_t2/bd_79_norm_t2.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_73
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_TP,
			iPS_Ctrl_TP_vs_iPS_BD_R_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
bd_73 = Reduce(intersect, x)		
length(bd_73)
write.table(bd_73, "bd_tbl_norm_t2/bd_73_norm_t2.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_16
x = list(iPS_Ctrl_TP_vs_iPS_BD_NR_TP, iPS_Ctrl_TP_vs_iPS_BD_TP, UnChanged_iPS_Ctrl_TP_vs_Primary_TP)
xInt = Reduce(intersect, x)
y = iPS_Ctrl_TP_vs_iPS_BD_R_TP
bd_16 = setdiff(xInt,y)
length(bd_16)
write.table(bd_16, "bd_tbl_norm_t2/bd_16_norm_t2.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_88
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_R_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
yInt = Reduce(union,y)
bd_88 = setdiff(xInt,yInt)
length(bd_88)
write.table(bd_88, "bd_tbl_norm_t2/bd_88_norm_t2.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_r_45
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
yInt = Reduce(union,y)
bd_r_45 = setdiff(xInt,yInt)
length(bd_r_45)
write.table(bd_r_45, "bd_tbl_norm_t2/bd_r_45_norm_t2.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_nr_60
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
yInt = Reduce(union,y)
bd_nr_60 = setdiff(xInt,yInt)
length(bd_nr_60)
write.table(bd_nr_60, "bd_tbl_norm_t2/bd_nr_60_norm_t2.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)



