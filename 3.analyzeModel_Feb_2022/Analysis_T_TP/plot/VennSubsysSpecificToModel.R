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

pdf("vennSubsysSpecificToModel.pdf")

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

Subsys_Primary_T <- iAstro_Primary_T$Var2
Subsys_iPS_Ctrl_T <- iAstro_iPS_Ctrl_T$Var2
Subsys_iPS_BD_T <- iAstro_iPS_BD_T$Var2
Subsys_iPS_BD_R_T <- iAstro_iPS_BD_R_T$Var2
Subsys_iPS_BD_NR_T <- iAstro_iPS_BD_NR_T$Var2

Subsys_Primary_TP <- iAstro_Primary_TP$Var2
Subsys_iPS_Ctrl_TP <- iAstro_iPS_Ctrl_TP$Var2
Subsys_iPS_BD_TP <- iAstro_iPS_BD_TP$Var2
Subsys_iPS_BD_R_TP <- iAstro_iPS_BD_R_TP$Var2
Subsys_iPS_BD_NR_TP <- iAstro_iPS_BD_NR_TP$Var2


################
# Subsys slice
################

#Subsys_Primary_TP
x = list(Subsys_iPS_Ctrl_TP, Subsys_iPS_BD_TP, Subsys_iPS_BD_R_TP, Subsys_iPS_BD_NR_TP)
xInt = Reduce(union, x)
y = Subsys_Primary_TP
Primary_TP = setdiff(y,xInt)
length(Primary_TP)
write.table(Primary_TP, "Subsys_Primary_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#Subsys_iPS_Ctrl_TP
x = list(Subsys_Primary_TP, Subsys_iPS_BD_TP, Subsys_iPS_BD_R_TP, Subsys_iPS_BD_NR_TP)
xInt = Reduce(union, x)
y = Subsys_iPS_Ctrl_TP
iPS_Ctrl_TP = setdiff(y,xInt)
length(iPS_Ctrl_TP)
write.table(iPS_Ctrl_TP, "Subsys_iPS_Ctrl_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#Subsys_iPS_BD_TP
x = list(Subsys_Primary_TP, Subsys_iPS_Ctrl_TP, Subsys_iPS_BD_R_TP, Subsys_iPS_BD_NR_TP)
xInt = Reduce(union, x)
y = Subsys_iPS_BD_TP
iPS_BD_TP = setdiff(y,xInt)
length(iPS_BD_TP)
write.table(iPS_BD_TP, "Subsys_iPS_BD_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#Subsys_iPS_BD_R_TP
x = list(Subsys_Primary_TP, Subsys_iPS_Ctrl_TP, Subsys_iPS_BD_TP, Subsys_iPS_BD_NR_TP)
xInt = Reduce(union, x)
y = Subsys_iPS_BD_R_TP
iPS_BD_R_TP = setdiff(y,xInt)
length(iPS_BD_R_TP)
write.table(iPS_BD_R_TP, "Subsys_iPS_BD_R_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#Subsys_iPS_BD_NR_TP
x = list(Subsys_Primary_TP, Subsys_iPS_Ctrl_TP, Subsys_iPS_BD_TP, Subsys_iPS_BD_R_TP)
xInt = Reduce(union, x)
y = Subsys_iPS_BD_NR_TP
iPS_BD_NR_TP = setdiff(y,xInt)
length(iPS_BD_NR_TP)
write.table(iPS_BD_NR_TP, "Subsys_iPS_BD_NR_TP.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

