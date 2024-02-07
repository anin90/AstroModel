library("pheatmap")
library("RColorBrewer")
library("ggplot2")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("ggpubr")
#library("ggpubr")
col.pal <- brewer.pal(9,"Spectral")

dat = read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/AstroModel/Databases_Tools/cobratoolbox/AstroModel_QC/2. sanityCheck/7_3_4/plot_sol/T3.csv", sep=",", header=T);
row.names(dat) <- dat$Objective
dat <- dat[,2:16]
dat_matrix <- data.matrix(dat)
dat_matrix[which(!is.finite(dat_matrix))] <- 3
dat_matrix = ifelse(dat_matrix>2,1,0)
dat_matrix = colSums(dat_matrix)

pdf("plots.pdf")
barplot(dat_matrix,las=2, cex.axis=1, col="blue")
dev.off()
