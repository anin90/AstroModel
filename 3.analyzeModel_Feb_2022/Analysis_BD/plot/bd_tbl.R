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

pdf("bd_tbl.pdf")

###########
# Load data
###########

bd_209_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_209_tbl.csv",
		header = T, sep = "\t")
bd_r_72_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_r_72_tbl.csv",
		header = T, sep = "\t")		
bd_nr_726_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_nr_726_tbl.csv",
		header = T, sep = "\t")

###########
# Plot data
###########
		
# Fig.1: venn rxns
	x <- list(
		BD_Lumped = bd_209_tbl$rxnList, 
		BD_Responder = bd_r_72_tbl$rxnList, 
		BD_NonResponder = bd_nr_726_tbl$rxnList)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlap in rxns disrupted between models")

# Fig.2: venn subsystems
	x <- list(
		BD_Lumped = bd_209_tbl$subSystem, 
		BD_Responder = bd_r_72_tbl$subSystem, 
		BD_NonResponder = bd_nr_726_tbl$subSystem)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlap in subSystems disrupted between models")				

# Fig.3: Number of rxns per subSystem (bd_lumped)
	iAstro_iPS_BD_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_TP.csv", header = T, sep = "\t")
	iAstro_iPS_BD_TP_subSystem = iAstro_iPS_BD_TP %>% count(Var2)
	bd_209_tbl_subSystem = bd_209_tbl %>% count(subSystem)
	comb = merge(bd_209_tbl_subSystem, iAstro_iPS_BD_TP_subSystem, by.x = "subSystem", by.y = "Var2")
	comb = comb[order(comb$n.x),]
	#3a
	dotchart(comb$n.x, labels = comb$subSystem,
			cex = 0.7, bg = "blue", xlab = "# of rxns", main = "# of disrupted rxns (n=209, bd_lumped)")
	#3b
	dotchart(comb$n.y, labels = comb$subSystem,
			cex = 0.7, bg = "grey", xlab = "# of rxns", main = "# of disrupted rxns (n=209, bd_lumped) vs all rxns")	
	points(comb$n.x, 1:nrow(comb), col = "blue", cex = 0.7)	
	#3c
	p.val = phyper(q = comb$n.x-1, m = comb$n.y, n = nrow(iAstro_iPS_BD_TP)-comb$n.y, k = sum(comb$n.x), lower.tail = FALSE, log.p = FALSE)
	p.val.fdr = p.adjust(p.val, method = "fdr")
	combStat = cbind(comb, p.val, p.val.fdr)
	combStat = combStat[order(-combStat$p.val.fdr),]
	combStat$subSystem <- factor(combStat$subSystem, levels = combStat$subSystem)
	ggplot(combStat) +
		geom_point(aes(y = subSystem, x = p.val.fdr),color='blue') +
		geom_vline(xintercept = 0.05, lty = 2) + theme_classic()+
		ggtitle("over-representation analysis, bd_lumped") +
		xlab("hypergeometric significance (fdr.adj.p.value)")
	#write fdr.significant subset to csv	
	bd_lumped_subsystem_fdr = subset(combStat, p.val.fdr <= 0.05)
	write.table(bd_lumped_subsystem_fdr, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_lumped_subsystem_fdr.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	bd_lumped_rxns_fdr = bd_209_tbl[bd_209_tbl$subSystem %in% bd_lumped_subsystem_fdr$subSystem,]
	write.table(bd_lumped_rxns_fdr, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_lumped_rxns_fdr.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
		
# Fig.4: Number of rxns per subSystem (bd_responder)
	iAstro_iPS_BD_R_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_R_TP.csv", header = T, sep = "\t")
	iAstro_iPS_BD_R_TP_subSystem = iAstro_iPS_BD_R_TP %>% count(Var2)
	bd_r_72_tbl_subSystem = bd_r_72_tbl %>% count(subSystem)
	comb = merge(bd_r_72_tbl_subSystem, iAstro_iPS_BD_R_TP_subSystem, by.x = "subSystem", by.y = "Var2")
	comb = comb[order(comb$n.x),]
	#4a
	dotchart(comb$n.x, labels = comb$subSystem,
			cex = 0.7, bg = "green", xlab = "# of rxns", main = "# of disrupted rxns (n=72, bd_responder)")
	#4b
	dotchart(comb$n.y, labels = comb$subSystem,
			cex = 0.7, bg = "grey", xlab = "# of rxns", main = "# of disrupted rxns (n=72, bd_responder) vs all rxns")	
	points(comb$n.x, 1:nrow(comb), col = "green", cex = 0.7)	
	#4c
	p.val = phyper(q = comb$n.x-1, m = comb$n.y, n = nrow(iAstro_iPS_BD_R_TP)-comb$n.y, k = sum(comb$n.x), lower.tail = FALSE, log.p = FALSE)
	p.val.fdr = p.adjust(p.val, method = "fdr")
	combStat = cbind(comb, p.val, p.val.fdr)
	combStat = combStat[order(-combStat$p.val.fdr),]
	combStat$subSystem <- factor(combStat$subSystem, levels = combStat$subSystem)
	ggplot(combStat) +
		geom_point(aes(y = subSystem, x = p.val.fdr),color='green') +
		geom_vline(xintercept = 0.05, lty = 2) + theme_classic()+
		ggtitle("over-representation analysis, bd_responder") +
		xlab("hypergeometric significance (fdr.adj.p.value)")
	#write fdr.significant subset to csv	
	bd_r_subsystem_fdr = subset(combStat, p.val.fdr <= 0.05)
	write.table(bd_r_subsystem_fdr, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_r_subsystem_fdr.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA) 
	bd_r_rxns_fdr = bd_r_72_tbl[bd_r_72_tbl$subSystem %in% bd_r_subsystem_fdr$subSystem,]
	write.table(bd_r_rxns_fdr, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_r_rxns_fdr.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

# Fig.5: Number of rxns per subSystem (bd_nonresponder)
	iAstro_iPS_BD_NR_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_NR_TP.csv", header = T, sep = "\t")
	iAstro_iPS_BD_NR_TP_subSystem = iAstro_iPS_BD_NR_TP %>% count(Var2)
	bd_nr_726_tbl_subSystem = bd_nr_726_tbl %>% count(subSystem)
	comb = merge(bd_nr_726_tbl_subSystem, iAstro_iPS_BD_NR_TP_subSystem, by.x = "subSystem", by.y = "Var2")
	comb = comb[order(comb$n.x),]
	#5a
	dotchart(comb$n.x, labels = comb$subSystem,
			cex = 0.7, bg = "red", xlab = "# of rxns", main = "# of disrupted rxns (n=726, bd_nonresponder)")	
	#5b
	dotchart(comb$n.y, labels = comb$subSystem,
			cex = 0.7, bg = "grey", xlab = "# of rxns", main = "# of disrupted rxns (n=726, bd_nonresponder) vs all rxns")	
	points(comb$n.x, 1:nrow(comb), col = "red", cex = 0.7)	
	#5c
	p.val = phyper(q = comb$n.x-1, m = comb$n.y, n = nrow(iAstro_iPS_BD_NR_TP)-comb$n.y, k = sum(comb$n.x), lower.tail = FALSE, log.p = FALSE)
	p.val.fdr = p.adjust(p.val, method = "fdr")
	combStat = cbind(comb, p.val, p.val.fdr)
	combStat = combStat[order(-combStat$p.val.fdr),]
	combStat$subSystem <- factor(combStat$subSystem, levels = combStat$subSystem)
	ggplot(combStat) +
		geom_point(aes(y = subSystem, x = p.val.fdr),color='red') +
		geom_vline(xintercept = 0.05, lty = 2) + theme_classic() +
		ggtitle("over-representation analysis, bd_nonresponder") +
		xlab("hypergeometric significance (fdr.adj.p.value)")	
	#write fdr.significant subset to csv	
	bd_nr_subsystem_fdr = subset(combStat, p.val.fdr <= 0.05)
	write.table(bd_nr_subsystem_fdr, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_nr_subsystem_fdr.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	bd_nr_rxns_fdr = bd_nr_726_tbl[bd_nr_726_tbl$subSystem %in% bd_nr_subsystem_fdr$subSystem,]
	write.table(bd_nr_rxns_fdr, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_nr_rxns_fdr.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

# Fig.6: Number of rxns per compartment (bd_lumped)
	#6a: all.209
	myTable <- table(bd_209_tbl$Var7, bd_209_tbl$Var9)
	compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "endoplasmic reticulum", "peroxisome", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") +
		scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
		labs(fill = "Direction of disruption", x ="compartment", y = "# of rxns") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
		theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
		scale_x_discrete(labels= compartment) + ggtitle("Disrupted rxns (all, n=209) per compartment (bd_lumped)")

	#6b: only.fdr.significant	
	myTable <- table(bd_lumped_rxns_fdr$Var7, bd_lumped_rxns_fdr$Var9)
	compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "endoplasmic reticulum", "peroxisome", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") +
		scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
		labs(fill = "Direction of disruption", x ="compartment", y = "# of rxns") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
		theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
		scale_x_discrete(labels= compartment) + ggtitle("Disrupted rxns (fdr.significant, n=90) per compartment (bd_lumped)")

# Fig.7: Number of rxns per compartment (bd_responder)
	#7a: all.72
	myTable <- table(bd_r_72_tbl$Var7, bd_r_72_tbl$Var9)
	compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "nucleus", "endoplasmic reticulum", 
					"peroxisome", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") +
		scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
		labs(fill = "Direction of disruption", x ="compartment", y = "# of rxns") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
		theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
		scale_x_discrete(labels= compartment) + ggtitle("Disrupted rxns (all, n=72) per compartment (bd_responder)")

	#7b: only.fdr.significant	
	myTable <- table(bd_r_rxns_fdr$Var7, bd_r_rxns_fdr$Var9)
	compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "nucleus", "endoplasmic reticulum", 
					"peroxisome", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") +
		scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
		labs(fill = "Direction of disruption", x ="compartment", y = "# of rxns") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
		theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
		scale_x_discrete(labels= compartment) + ggtitle("Disrupted rxns (fdr.significant, n=36) per compartment (bd_responder)")

# Fig.8: Number of rxns per compartment (bd_nonresponder)
	#8a: all.726
	myTable <- table(bd_nr_726_tbl$Var7, bd_nr_726_tbl$Var9)
	compartment <- c("cytoplasm", "extracellular space", "Golgi apparatus", "lysosome", "mitochondrion", "nucleus", 
					"endoplasmic reticulum", "peroxisome", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") +
		scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
		labs(fill = "Direction of disruption", x ="compartment", y = "# of rxns") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
		theme(axis.text.x = element_text(size=12, angle = 60, hjust=1))+ 
		scale_x_discrete(labels= compartment) + ggtitle("Disrupted rxns (all, n=726) per compartment (bd_nonresponder)")

	#8b: only.fdr.significant	
	myTable <- table(bd_nr_rxns_fdr$Var7, bd_nr_rxns_fdr$Var9)
	compartment <- c("cytoplasm", "extracellular space", "Golgi apparatus", "lysosome", "mitochondrion", "nucleus", 
					"endoplasmic reticulum", "peroxisome", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") +
		scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
		labs(fill = "Direction of disruption", x ="compartment", y = "# of 'rxns") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
		theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
		scale_x_discrete(labels= compartment) + ggtitle("Disrupted rxns (fdr.significant, n=263) per compartment (bd_nonresponder)")
		
# Fig.9: venn rxns (fdr.significant)
	x <- list(
		BD_Lumped = bd_lumped_rxns_fdr$rxnList, 
		BD_Responder = bd_r_rxns_fdr$rxnList, 
		BD_NonResponder = bd_nr_rxns_fdr$rxnList)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlp in rxns (fdr.significant) disrupted between models")

# Fig.10: venn subSystems (fdr.significant)
	x <- list(
		BD_Lumped = bd_lumped_rxns_fdr$subSystem, 
		BD_Responder = bd_r_rxns_fdr$subSystem, 
		BD_NonResponder = bd_nr_rxns_fdr$subSystem)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlp in subSystems (fdr.significant) disrupted between models")		
