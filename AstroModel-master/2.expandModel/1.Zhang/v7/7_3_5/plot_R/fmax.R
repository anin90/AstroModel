library("pheatmap")
library("RColorBrewer")
dat = read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/AstroModel/Databases_Tools/cobratoolbox/AstroModel_QC/2. sanityCheck/7_3_4/plot_sol/T3.csv", sep=",", header=T);
row.names(dat) <- dat$Objective
dat <- dat[,2:16]
dat_matrix <- data.matrix(dat)
dat_matrix[which(!is.finite(dat_matrix))] <- 2.001
dat_matrix = ifelse(dat_matrix>2,1,0)

#col.pal <- brewer.pal(9,"RdYlBu")
col.pal <- brewer.pal(9,"Blues")
hm.parameters <- list(dat_matrix,color = col.pal, 
fontsize_col = 7, fontsize_row = 5, scale = "none",
cluster_rows = TRUE, cluster_cols = TRUE,border_color = "grey60",cellwidth = 15)

#do.call("pheatmap", hm.parameters)
#dev.print(pdf, 'T1.pdf')

pdf(file="T3.pdf")
do.call("pheatmap", hm.parameters)
dev.off()
