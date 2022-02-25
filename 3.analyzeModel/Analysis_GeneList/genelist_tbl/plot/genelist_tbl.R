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

pdf("genelist_tbl.pdf", width=7, height=2)

###########
# Load data
###########

CKMT2_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/CKMT2_Primary_TP_tbl.csv", header = T, sep = "\t")
INPP5A_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/INPP5A_Primary_TP_tbl.csv", header = T, sep = "\t")
NMNAT1_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/NMNAT1_Primary_TP_tbl.csv", header = T, sep = "\t")
SLC5A8_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/SLC5A8_Primary_TP_tbl.csv", header = T, sep = "\t")
PLA2G12A_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/PLA2G12A_Primary_TP_tbl.csv", header = T, sep = "\t")
SLC29A2_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/SLC29A2_Primary_TP_tbl.csv", header = T, sep = "\t")
COASY_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/COASY_Primary_TP_tbl.csv", header = T, sep = "\t")
SYNJ2_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/SYNJ2_Primary_TP_tbl.csv", header = T, sep = "\t")
PYGL_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/PYGL_Primary_TP_tbl.csv", header = T, sep = "\t")
SLC22A9_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/SLC22A9_Primary_TP_tbl.csv", header = T, sep = "\t")
PLCB4_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/PLCB4_Primary_TP_tbl.csv", header = T, sep = "\t")
PLA2G6_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/PLA2G6_Primary_TP_tbl.csv", header = T, sep = "\t")
GBA_Primary <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/GBA_Primary_TP_tbl.csv", header = T, sep = "\t")


CKMT2_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/CKMT2_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
INPP5A_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/INPP5A_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
NMNAT1_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/NMNAT1_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
SLC5A8_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/SLC5A8_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
PLA2G12A_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/PLA2G12A_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
SLC29A2_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/SLC29A2_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
COASY_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/COASY_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
SYNJ2_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/SYNJ2_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
PYGL_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/PYGL_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
SLC22A9_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/SLC22A9_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
PLCB4_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/PLCB4_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
PLA2G6_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/PLA2G6_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")
GBA_iPS_Ctrl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_GeneList/genelist_tbl/GBA_iPS_Ctrl_TP_tbl.csv", header = T, sep = "\t")

######################
# Plot data- for loop
######################

#~ models<- list(CKMT2_Primary, CKMT2_iPS_Ctrl, INPP5A_Primary, INPP5A_iPS_Ctrl, NMNAT1_Primary, NMNAT1_iPS_Ctrl)
#~ modelnames <- c('CKMT2_Primary', 'CKMT2_iPS_Ctrl', 'INPP5A_Primary', 'INPP5A_iPS_Ctrl', 'NMNAT1_Primary', 'NMNAT1_iPS_Ctrl')

#~ for (i in models) {
#~ 		dat = i
#~ 		myTable = table(dat$PerturbationClass, dat$subSystem)
#~ 		print(ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
#~ 			geom_bar(stat="identity") + coord_flip() +
#~ 			scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
#~ 			labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
#~ 			geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)))+
#~ 			ggtitle(paste0(unlist(i)))
#~ 		}

###########
# Plot data
###########

# Fig.1a: CKMT2
	dat = CKMT2_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2019 - CKMT2_Primary")
		
	# Fig.1b:
	dat = CKMT2_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2019 - CKMT2_iPS_Ctrl")

# Fig.2a: INPP5A
	dat = INPP5A_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2019 - INPP5A_Primary")
		
	# Fig.2b:
	dat = INPP5A_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2019 - INPP5A_iPS_Ctrl")

# Fig.3a: NMNAT1
	dat = NMNAT1_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - NMNAT1_Primary")
		
	# Fig.3b:
	dat = NMNAT1_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - NMNAT1_iPS_Ctrl")

# Fig.4a: SLC5A8
	dat = SLC5A8_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - SLC5A8_Primary")
		
	# Fig.4b:
	dat = SLC5A8_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - SLC5A8_iPS_Ctrl")

# Fig.5a: PLA2G12A
	dat = PLA2G12A_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - PLA2G12A_Primary")
		
	# Fig.5b:
	dat = PLA2G12A_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - PLA2G12A_iPS_Ctrl")

# Fig.6a: SLC29A2
	dat = SLC29A2_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - SLC29A2_Primary")
		
	# Fig.6b:
	dat = SLC29A2_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - SLC29A2_iPS_Ctrl")

# Fig.7a: COASY
	dat = COASY_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - COASY_Primary")
		
	# Fig.7b:
	dat = COASY_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - COASY_iPS_Ctrl")

# Fig.8a: SYNJ2
	dat = SYNJ2_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - SYNJ2_Primary")
		
	# Fig.8b:
	dat = SYNJ2_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - SYNJ2_iPS_Ctrl")

# Fig.9a: PYGL
	dat = PYGL_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - PYGL_Primary")
		
	# Fig.9b:
	dat = PYGL_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - PYGL_iPS_Ctrl")

# Fig.10a: SLC22A9
	dat = SLC22A9_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - SLC22A9_Primary")
		
	# Fig.10b:
	dat = SLC22A9_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Suhas 2021 - SLC22A9_iPS_Ctrl")

# Fig.11a: PLCB4
	dat = PLCB4_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Syama 2018 - PLCB4_Primary")
		
	# Fig.11b:
	dat = PLCB4_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Syama 2018 - PLCB4_iPS_Ctrl")

# Fig.12a: PLA2G6
	dat = PLA2G6_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Booth 2017 - PLA2G6_Primary")
		
	# Fig.12b:
	dat = PLA2G6_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Booth 2017 - PLA2G6_iPS_Ctrl")

# Fig.13a: GBA
	dat = GBA_Primary
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Booth 2017 - GBA_Primary")
		
	# Fig.13b:
	dat = GBA_iPS_Ctrl
	head(dat)
	myTable = table(dat$PerturbationClass, dat$subSystem)
	head(myTable)
	ggplot(as.data.frame(myTable), aes(x = reorder(Var2, -Freq), Freq, fill=Var1, label = Freq))+
		geom_bar(stat="identity") + coord_flip() +
		scale_fill_manual(values = c("khaki", "pink1")) + theme_classic() +
		labs(fill = "PerturbationClass", x ="subsystem", y = "# of disrupted reactions") +
		geom_text(data=subset(as.data.frame(myTable),Freq != 0), position = position_stack(vjust = 0.5)) +
		ggtitle("Booth 2017 - GBA_iPS_Ctrl")
