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

pdf("bd_tbl_significant_norm_t2/plotFinalTable_norm_t2.pdf", width=7, height=3)

###########
# Load data
###########

# load subsystems
bd_lumped_subsystem_fdr <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant_norm_t2/bd_lumped_subsystem_fdr_norm_t2.csv",
		header = T, sep = "\t")
bd_r_subsystem_fdr <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant_norm_t2/bd_r_subsystem_fdr_norm_t2.csv",
		header = T, sep = "\t")		
bd_nr_subsystem_fdr <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant_norm_t2/bd_nr_subsystem_fdr_norm_t2.csv",
		header = T, sep = "\t")

# load rxns
bd_lumped_rxns_fdr <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant_norm_t2/bd_lumped_rxns_fdr_norm_t2.csv",
		header = T, sep = "\t")
bd_r_rxns_fdr <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant_norm_t2/bd_r_rxns_fdr_norm_t2.csv",
		header = T, sep = "\t")		
bd_nr_rxns_fdr <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant_norm_t2/bd_nr_rxns_fdr_norm_t2.csv",
		header = T, sep = "\t")

###########
# Plot data
###########

# Fig.1a: Number of rxns per subSystem post fdr (BD_Lumped)
	mdat = bd_lumped_subsystem_fdr
	ggplot(data=mdat, aes(x = reorder(subSystem, -n.x), n.x))+ theme_classic() +
		theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
		geom_bar(stat="identity",aes(fill = p.val.fdr)) + coord_flip() + 
		labs(fill = "p.val.fdr (<=0.05)", x ="subsystem", y = "# of disrupted reactions") +
		ggtitle("BD_Lumped")

	# Fig.1b:
	myTable <- table(bd_lumped_rxns_fdr$Flux, bd_lumped_rxns_fdr$subSystem)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+ 
		geom_bar(stat="identity") + coord_flip() + 
		scale_fill_manual(values = c("palegreen3", "orangered1")) + theme_classic() + 
		labs(fill = "Flux", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("BD_Lumped")	
		
# Fig.2: Number of rxns per subSystem post fdr (BD_Responder)
	mdat = bd_r_subsystem_fdr
	ggplot(data=mdat, aes(x = reorder(subSystem, -n.x), n.x))+ theme_classic() +
		theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
		geom_bar(stat="identity",aes(fill = p.val.fdr), width = 0.5) + coord_flip() + 
		labs(fill = "p.val.fdr (<=0.05)", x ="subsystem", y = "# of disrupted reactions") +
		ggtitle("BD_Responder")
		
	# Fig.2b:
	myTable <- table(bd_r_rxns_fdr$Flux, bd_r_rxns_fdr$subSystem)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+ 
		geom_bar(stat="identity", width = 0.5) + coord_flip() + 
		scale_fill_manual(values = c("palegreen3", "orangered1")) + theme_classic() + 
		labs(fill = "Flux", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("BD_Responder")	
		
# Fig.3: Number of rxns per subSystem post fdr (BD_NonResponder)
	mdat = bd_nr_subsystem_fdr
	ggplot(data=mdat, aes(x = reorder(subSystem, -n.x), n.x))+ theme_classic() +
		theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))+
		geom_bar(stat="identity",aes(fill = p.val.fdr)) + coord_flip() + 
		labs(fill = "p.val.fdr (<=0.05)", x ="subsystem", y = "# of disrupted reactions") +
		ggtitle("BD_NonResponder") 

	# Fig.3b:
	myTable <- table(bd_nr_rxns_fdr$Flux, bd_nr_rxns_fdr$subSystem)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+ 
		geom_bar(stat="identity") + coord_flip() + 
		scale_fill_manual(values = c("palegreen3", "orangered1")) + theme_classic() + 
		labs(fill = "Flux", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("BD_NonResponder")						
				
# Fig.4: Number of rxns per compartment post fdr (BD_Lumped)
	myTable <- table(bd_lumped_rxns_fdr$Flux, bd_lumped_rxns_fdr$Localization)
	head(myTable)
	compartment <- c("cytoplasm", "lysosome", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") + coord_flip() + 
		scale_fill_manual(values = c("palegreen3", "orangered1")) + theme_classic() + 
		labs(fill = "Flux", x ="compartment", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		scale_x_discrete(labels= compartment) + 
		ggtitle("BD_Lumped")	
		
# Fig.5: Number of rxns per compartment post fdr (BD_Responder)
	myTable <- table(bd_r_rxns_fdr$Flux, bd_r_rxns_fdr$Localization)
	head(myTable)	
	compartment <- c("cytoplasm", "Golgi apparatus", "lysosome", "endoplasmic reticulum", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") + coord_flip() + 
		scale_fill_manual(values = c("palegreen3", "orangered1")) + theme_classic() + 
		labs(fill = "Flux", x ="compartment", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		scale_x_discrete(labels= compartment) + 
		ggtitle("BD_Responder")		

# Fig.6: Number of rxns per compartment post fdr (BD_NonResponder)
	myTable <- table(bd_nr_rxns_fdr$Flux, bd_nr_rxns_fdr$Localization)
	head(myTable)	
	compartment <- c("cytoplasm", "extracellular space", "mitochondrion", "peroxisome", "intercompartmental")
	ggplot(as.data.frame(myTable), aes(Var2, Freq, fill=Var1, label = Freq)) + 
		geom_bar(stat="identity") + coord_flip() + 
		scale_fill_manual(values = c("palegreen3", "orangered1")) + theme_classic() + 
		labs(fill = "Flux", x ="compartment", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		scale_x_discrete(labels= compartment) + 
		ggtitle("BD_NonResponder")		

