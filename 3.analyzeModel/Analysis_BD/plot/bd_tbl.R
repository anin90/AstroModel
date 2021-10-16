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
library(plotrix)

pdf("bd_tbl.pdf")

###########
# Load data
###########

bd_212_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_212_tbl.csv",
		header = T, sep = "\t")
bd_r_92_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_r_92_tbl.csv",
		header = T, sep = "\t")		
bd_nr_670_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_nr_670_tbl.csv",
		header = T, sep = "\t")
		
# Fig.1: venn rxns
x <- list(
  BD_Lumped = bd_212_tbl$rxnList, 
  BD_Responder = bd_r_92_tbl$rxnList, 
  BD_NonResponder = bd_nr_670_tbl$rxnList)
ggvenn(x, fill_color = c("blue", "green", "red"))

# Fig.2: venn subsystems
x <- list(
  BD_Lumped = bd_212_tbl$subSystem, 
  BD_Responder = bd_r_92_tbl$subSystem, 
  BD_NonResponder = bd_nr_670_tbl$subSystem)
ggvenn(x, fill_color = c("blue", "green", "red"))				


# Fig.3: Number of rxns per subSystem (bd_lumped)
iAstro_iPS_BD_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_TP_subSystem = iAstro_iPS_BD_TP %>% count(Var2)
bd_212_tbl_subSystem = bd_212_tbl %>% count(subSystem)
comb = merge(bd_212_tbl_subSystem, iAstro_iPS_BD_TP_subSystem, by.x = "subSystem", by.y = "Var2")
comb = comb[order(comb$n.x),]
#3a
dotchart(comb$n.x, labels = comb$subSystem,
			cex = 0.7, bg = "blue", xlab = "# of disrupted rxns (n=212, bd_lumped) per subSystem")	
#3b
dotchart(comb$n.y, labels = comb$subSystem,
			cex = 0.7, bg = "grey", xlab = "# of disrupted rxns (n=212, bd_lumped) per subSystem")	
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
    ggtitle("BD_Lumped metabolic model") +
  xlab("hypergeometric significance (fdr.adj.p.value)") + ylab("Subsystems disrupted")


# Fig.4: Number of rxns per subSystem (bd_responder)
iAstro_iPS_BD_R_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_R_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_R_TP_subSystem = iAstro_iPS_BD_R_TP %>% count(Var2)
bd_r_92_tbl_subSystem = bd_r_92_tbl %>% count(subSystem)
comb = merge(bd_r_92_tbl_subSystem, iAstro_iPS_BD_R_TP_subSystem, by.x = "subSystem", by.y = "Var2")
comb = comb[order(comb$n.x),]
#4a
dotchart(comb$n.x, labels = comb$subSystem,
			cex = 0.7, bg = "green", xlab = "# of disrupted rxns (n=92, bd_responder) per subSystem")	
#4b
dotchart(comb$n.y, labels = comb$subSystem,
			cex = 0.7, bg = "grey", xlab = "# of disrupted rxns (n=92, bd_responder) per subSystem")	
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
    ggtitle("BD_Responder metabolic model") +
  xlab("hypergeometric significance (fdr.adj.p.value)") + ylab("Subsystems disrupted")


# Fig.5: Number of rxns per subSystem (bd_nonresponder)
iAstro_iPS_BD_NR_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_NR_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_NR_TP_subSystem = iAstro_iPS_BD_NR_TP %>% count(Var2)
bd_nr_670_tbl_subSystem = bd_nr_670_tbl %>% count(subSystem)
comb = merge(bd_nr_670_tbl_subSystem, iAstro_iPS_BD_NR_TP_subSystem, by.x = "subSystem", by.y = "Var2")
comb = comb[order(comb$n.x),]
#5a
dotchart(comb$n.x, labels = comb$subSystem,
			cex = 0.7, bg = "red", xlab = "# of disrupted rxns (n=670, bd_nonresponder) per subSystem")	
#5b
dotchart(comb$n.y, labels = comb$subSystem,
			cex = 0.7, bg = "grey", xlab = "# of disrupted rxns (n=670, bd_nonresponder) per subSystem")	
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
    ggtitle("BD_Non-responder metabolic model") +
  xlab("hypergeometric significance (fdr.adj.p.value)") + ylab("Subsystems disrupted")	



# Fig.6: Number of rxns per compartment (bd_lumped)
myTable <- table(bd_212_tbl$Var7, bd_212_tbl$Var8)
compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "endoplasmic reticulum", "intercompartmental")

ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
   labs(fill = "Direction of disruption", x ="Compartment", y = "# of disrupted rxns (n=212), bd_lumped)") +
	geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
	theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
	scale_x_discrete(labels= compartment)

# Fig.7: Number of rxns per compartment (bd_responder)
myTable <- table(bd_r_92_tbl$Var7, bd_r_92_tbl$Var8)
compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "nucleus", "endoplasmic reticulum", 
				"peroxisome", "intercompartmental")

ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
   labs(fill = "Direction of disruption", x ="Compartment", y = "# of disrupted rxns (n=92), bd_responder)") +
	geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
	theme(axis.text.x = element_text(size=12, angle = 60, hjust=1)) + 
	scale_x_discrete(labels= compartment)
  
# Fig.8: Number of rxns per compartment (bd_nonresponder)
myTable <- table(bd_nr_670_tbl$Var7, bd_nr_670_tbl$Var8)
compartment <- c("cytoplasm", "extracellular space", "Golgi apparatus", "lysosome", "mitochondrion", "nucleus", 
				"endoplasmic reticulum", "peroxisome", "intercompartmental")

ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("green", "orange")) + theme_classic() + 
   labs(fill = "Direction of disruption", x ="Compartment", y = "# of disrupted rxns (n=670), bd_nonresponder)") +
	geom_text(data=subset(as.data.frame(myTable),Freq != 0), size = 3, position = position_stack(vjust = 0.5)) +
	theme(axis.text.x = element_text(size=12, angle = 60, hjust=1))+ 
	scale_x_discrete(labels= compartment)














