# http://www.sthda.com/english/wiki/dot-charts-r-base-graphs
# https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
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
library(plotrix)
library(pheatmap)
library(magrittr)
library(venn)

pdf("st_tbl_abs/st_tbl_abs.pdf")

###########
# Load data
###########

ST <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_abs/st_rxns.csv",
		header = T, sep = "\t")

###########
# Plot data
###########		

# Fig.1: Number of rxns per subSystem (st_rxns)
	iPS_ST_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_ST_TP_abs.csv", header = T, sep = "\t")
	iPS_ST_TP_subSystem = iPS_ST_TP %>% count(SubSystem)
	ST_subSystem = ST %>% count(subSystem)
	comb = merge(ST_subSystem, iPS_ST_TP_subSystem, by.x = "subSystem", by.y = "SubSystem")
	comb = comb[order(comb$n.x),]
	#1a
	dotchart(comb$n.x, labels = comb$subSystem,
			cex = 0.7, bg = "blue", xlab = "# of rxns", main = "# of disrupted rxns (n=120, ST)")
	#1b
	dotchart(comb$n.y, labels = comb$subSystem,
			cex = 0.7, bg = "grey", xlab = "# of rxns", main = "# of disrupted rxns (n=120, ST) vs all rxns")	
	points(comb$n.x, 1:nrow(comb), col = "blue", cex = 0.7)	
	#1c
	p.val = phyper(q = comb$n.x-1, m = comb$n.y, n = nrow(iPS_ST_TP)-comb$n.y, k = sum(comb$n.x), lower.tail = FALSE, log.p = FALSE)
	p.val.fdr = p.adjust(p.val, method = "fdr")
	combStat = cbind(comb, p.val, p.val.fdr)
	combStat = combStat[order(-combStat$p.val.fdr),]
	combStat$subSystem <- factor(combStat$subSystem, levels = combStat$subSystem)
	ggplot(combStat) +
		geom_point(aes(y = subSystem, x = p.val.fdr),color='blue') +
		geom_vline(xintercept = 0.05, lty = 2) + theme_classic()+
		ggtitle("over-representation analysis, ST") +
		xlab("hypergeometric significance (fdr.adj.p.value)")
	#write fdr.significant subset to csv	
	ST_subsystem_fdr = subset(combStat, p.val.fdr <= 0.05)
	write.table(ST_subsystem_fdr, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_significant_abs/ST_subsystem_fdr_abs.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	ST_rxns_fdr = ST[ST$subSystem %in% ST_subsystem_fdr$subSystem,]
	write.table(ST_rxns_fdr, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_significant_abs/ST_rxns_fdr_abs.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)

# Fig.2: Number of rxns per compartment (st_rxns)
	#2a: all.120
	myTable <- table(ST$Flux, ST$Localization)
	compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "endoplasmic reticulum", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") +
		scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
		labs(fill = "Flux", x ="compartment", y = "# of rxns") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
		theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
		scale_x_discrete(labels= compartment) + ggtitle("Disrupted rxns (all, n=120) per compartment (ST)")

	#2b: only.fdr.significant	
	myTable <- table(ST_rxns_fdr$Flux, ST_rxns_fdr$Localization)	
	compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "endoplasmic reticulum", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") +
		scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
		labs(fill = "Flux", x ="compartment", y = "# of rxns") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
		theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
		scale_x_discrete(labels= compartment) + ggtitle("Disrupted rxns (fdr.significant, n=17) per compartment (ST)")
