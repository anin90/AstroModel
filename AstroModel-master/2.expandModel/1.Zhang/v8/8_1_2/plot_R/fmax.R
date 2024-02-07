library("pheatmap")
library("RColorBrewer")
dat = read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/AstroModel/Analysis/cobratoolbox/AstroModel/2_expandModel/8_1_2/plot_R/T11_plot.csv", sep=",", header=T);
row.names(dat) <- dat$Li_action
dat <- dat[,2:4]
dat_matrix <- data.matrix(dat)

#col.pal <- brewer.pal(9,"RdYlBu")
col.pal <- brewer.pal(9,"Reds")
hm.parameters <- list(dat_matrix,color = col.pal,cluster_rows = FALSE, cluster_cols = FALSE, 
border_color = "grey60",cellwidth = 18,cellheight = 18, scale = "none",legend = TRUE,
legend_breaks = c(0,1.5), legend_labels = c("1 (FSr_WT)",">1.5 (FSr_Li+)"))

 
pdf(file="T3.pdf")
do.call("pheatmap", hm.parameters)
dev.off()
