library("pheatmap")
library("RColorBrewer")
dat = read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/AstroModel/Analysis/cobratoolbox/AstroModel/2_expandModel/8_1_2/plot_R/v8_1_2_Tests.csv", sep=",", header=T);
row.names(dat) <- dat$Objective

boxplot(max_cTv~Objective_class,data=dat,main="Model function tests",
xlab="Metabolic objectives",ylab="Flux (mmol/gDW/h)",col = c("lightgreen","khaki1","salmon1"),
border="brown",las=1,outline=FALSE, boxwex=0.5, cex.main=1.5, cex.lab=1.2, cex.axis=1.15)
