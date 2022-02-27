# http://www.sthda.com/english/wiki/dot-charts-r-base-graphs
# https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
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
library(UpSetR)

pdf("vennRxnSpecificToModel.pdf")

###########
# Load data
###########

iAstro_Primary_T <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_Primary_T.csv", header = T, sep = "\t")
iAstro_iPS_Ctrl_T <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_Ctrl_T.csv", header = T, sep = "\t")
iAstro_iPS_BD_T <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_T.csv", header = T, sep = "\t")
iAstro_iPS_BD_R_T <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_R_T.csv", header = T, sep = "\t")
iAstro_iPS_BD_NR_T <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_NR_T.csv", header = T, sep = "\t")

iAstro_Primary_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_Primary_TP.csv", header = T, sep = "\t")
iAstro_iPS_Ctrl_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_Ctrl_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_R_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_R_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_NR_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_NR_TP.csv", header = T, sep = "\t")

Rxns_Primary_T <- iAstro_Primary_T$Var1
Rxns_iPS_Ctrl_T <- iAstro_iPS_Ctrl_T$Var1
Rxns_iPS_BD_T <- iAstro_iPS_BD_T$Var1
Rxns_iPS_BD_R_T <- iAstro_iPS_BD_R_T$Var1
Rxns_iPS_BD_NR_T <- iAstro_iPS_BD_NR_T$Var1

Rxns_Primary_TP <- iAstro_Primary_TP$Var1
Rxns_iPS_Ctrl_TP <- iAstro_iPS_Ctrl_TP$Var1
Rxns_iPS_BD_TP <- iAstro_iPS_BD_TP$Var1
Rxns_iPS_BD_R_TP <- iAstro_iPS_BD_R_TP$Var1
Rxns_iPS_BD_NR_TP <- iAstro_iPS_BD_NR_TP$Var1


###################################
# rnxs to matrix - wrangle
###################################
mat = lst(Rxns_Primary_TP, Rxns_iPS_Ctrl_TP, Rxns_iPS_BD_TP, Rxns_iPS_BD_R_TP, Rxns_iPS_BD_NR_TP,
			Rxns_Primary_T, Rxns_iPS_Ctrl_T, Rxns_iPS_BD_T, Rxns_iPS_BD_R_T, Rxns_iPS_BD_NR_T) %>% 			
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)
  
mat = t(mat)
write.table(mat, "Rxn_Table.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/plot/Rxn_Table.csv",
		header = T, sep = "\t")

##################################
# rnxs to matrix - UpSetR (TP)
##################################

p = upset(mat,
	order.by = c("degree","freq"),
	matrix.color="black", sets = c("Rxns_Primary_TP", "Rxns_iPS_Ctrl_TP", "Rxns_iPS_BD_TP", "Rxns_iPS_BD_R_TP", "Rxns_iPS_BD_NR_TP"), 
	keep.order = TRUE, sets.bar.color=c("dark grey","light grey","blue","green","red"))
p

################
# rnxs slice
################

#Rxns_Primary_TP
x = list(Rxns_iPS_Ctrl_TP, Rxns_iPS_BD_TP, Rxns_iPS_BD_R_TP, Rxns_iPS_BD_NR_TP)
xInt = Reduce(union, x)
y = Rxns_Primary_TP
Primary_TP = setdiff(y,xInt)
length(Primary_TP)
x <- subset(iAstro_Primary_TP, iAstro_Primary_TP$Var1 %in% Primary_TP)
write.table(x, "Rxns_Primary_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#Rxns_iPS_Ctrl_TP
x = list(Rxns_Primary_TP, Rxns_iPS_BD_TP, Rxns_iPS_BD_R_TP, Rxns_iPS_BD_NR_TP)
xInt = Reduce(union, x)
y = Rxns_iPS_Ctrl_TP
iPS_Ctrl_TP = setdiff(y,xInt)
length(iPS_Ctrl_TP)
x <- subset(iAstro_iPS_Ctrl_TP, iAstro_iPS_Ctrl_TP$Var1 %in% iPS_Ctrl_TP)
write.table(x, "Rxns_iPS_Ctrl_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#Rxns_iPS_BD_TP
x = list(Rxns_Primary_TP, Rxns_iPS_Ctrl_TP, Rxns_iPS_BD_R_TP, Rxns_iPS_BD_NR_TP)
xInt = Reduce(union, x)
y = Rxns_iPS_BD_TP
iPS_BD_TP = setdiff(y,xInt)
length(iPS_BD_TP)
x <- subset(iAstro_iPS_BD_TP, iAstro_iPS_BD_TP$Var1 %in% iPS_BD_TP)
write.table(x, "Rxns_iPS_BD_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#Rxns_iPS_BD_R_TP
x = list(Rxns_Primary_TP, Rxns_iPS_Ctrl_TP, Rxns_iPS_BD_TP, Rxns_iPS_BD_NR_TP)
xInt = Reduce(union, x)
y = Rxns_iPS_BD_R_TP
iPS_BD_R_TP = setdiff(y,xInt)
length(iPS_BD_R_TP)
x <- subset(iAstro_iPS_BD_R_TP, iAstro_iPS_BD_R_TP$Var1 %in% iPS_BD_R_TP)
write.table(x, "Rxns_iPS_BD_R_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#Rxns_iPS_BD_NR_TP
x = list(Rxns_Primary_TP, Rxns_iPS_Ctrl_TP, Rxns_iPS_BD_TP, Rxns_iPS_BD_R_TP)
xInt = Reduce(union, x)
y = Rxns_iPS_BD_NR_TP
iPS_BD_NR_TP = setdiff(y,xInt)
length(iPS_BD_NR_TP)
x <- subset(iAstro_iPS_BD_NR_TP, iAstro_iPS_BD_NR_TP$Var1 %in% iPS_BD_NR_TP)
write.table(x, "Rxns_iPS_BD_NR_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)



