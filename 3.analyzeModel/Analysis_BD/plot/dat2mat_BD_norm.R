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

pdf("FSR_UpSet_BD_norm.pdf")

##########################################
# rnxs to matrix: import output
##########################################

iPS_Ctrl_TP_vs_iPS_BD_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_norm_vs_iAstro_iPS_BD_TP_norm.csv",header = T, sep = "\t");
iPS_Ctrl_TP_vs_iPS_BD_R_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_norm_vs_iAstro_iPS_BD_R_TP_norm.csv",header = T, sep = "\t");
iPS_Ctrl_TP_vs_iPS_BD_NR_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_norm_vs_iAstro_iPS_BD_NR_TP_norm.csv",header = T, sep = "\t");
UnChanged_iPS_Ctrl_TP_vs_Primary_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Ctrl/UnChanged_iAstro_iPS_Ctrl_TP_norm_vs_iAstro_Primary_TP_norm.csv",header = T, sep = "\t");

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

write.table(mat, "FSR_Table_BD_norm.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/FSR_Table_BD_norm.csv",
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

svg(filename="FSR_UpSet_BD_TP_norm.svg")
p 
dev.off()

################
# rnxs slice
################

#bd_30
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
xInt = Reduce(intersect, x)
y = iPS_Ctrl_TP_vs_iPS_BD_NR_TP
bd_30 = setdiff(xInt,y)
length(bd_30)
write.table(bd_30, "bd_30_norm.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_38
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_TP,
			iPS_Ctrl_TP_vs_iPS_BD_R_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
bd_38 = Reduce(intersect, x)		
length(bd_38)
write.table(bd_38, "bd_38_norm.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_36
x = list(iPS_Ctrl_TP_vs_iPS_BD_NR_TP, iPS_Ctrl_TP_vs_iPS_BD_TP, UnChanged_iPS_Ctrl_TP_vs_Primary_TP)
xInt = Reduce(intersect, x)
y = iPS_Ctrl_TP_vs_iPS_BD_R_TP
bd_36 = setdiff(xInt,y)
length(bd_36)
write.table(bd_36, "bd_36_norm.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_95
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_R_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
yInt = Reduce(union,y)
bd_95 = setdiff(xInt,yInt)
length(bd_95)
write.table(bd_95, "bd_95_norm.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_r_155
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
yInt = Reduce(union,y)
bd_r_155 = setdiff(xInt,yInt)
length(bd_r_155)
write.table(bd_r_155, "bd_r_155_norm.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#bd_nr_148
x = list(UnChanged_iPS_Ctrl_TP_vs_Primary_TP, iPS_Ctrl_TP_vs_iPS_BD_NR_TP)
xInt = Reduce(intersect,x)
y = list(iPS_Ctrl_TP_vs_iPS_BD_TP, iPS_Ctrl_TP_vs_iPS_BD_R_TP)
yInt = Reduce(union,y)
bd_nr_148 = setdiff(xInt,yInt)
length(bd_nr_148)
write.table(bd_nr_148, "bd_nr_148_norm.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)
