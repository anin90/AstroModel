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
library(scales)
library(plyr)
pdf("Test4HumanAstro.pdf", width=3.5, height=3.5)

###########
# Load data
###########

Test4HumanAstro <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Test4HumanAstro_TP.csv",
		header = T, sep = "\t")
#~ Test4HumanAstro_T <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Test4HumanAstro_T.csv",
#~ 		header = T, sep = "\t")

# PI5P4K
	ss.data = Test4HumanAstro[Test4HumanAstro$TestSolutionName == 'PI5P4K',]						
	mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
            measure.vars=c("Primary", "iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
	mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))     
    mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl',]
	FluxCtrl = mdat_ctrl$value               

	ggplot(mdat, aes(x = variable, y = value)) +
		geom_col(aes(fill = TestSolutionName), width = 0.7, colour = "black") + theme_classic()+
		labs(fill = "Reaction ID             ", x ="modelOri", y = "Flux (mmol/gDw/h)")+
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		scale_fill_brewer(palette = "Blues") +
		geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","red"), size=1)

# MI1PP
	ss.data = Test4HumanAstro[Test4HumanAstro$TestSolutionName == 'MI1PP',]						
	mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
            measure.vars=c("Primary", "iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
	mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))     
    mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl',]
	FluxCtrl = mdat_ctrl$value               

	ggplot(mdat, aes(x = variable, y = value)) +
		geom_col(aes(fill = TestSolutionName), width = 0.7, colour = "black") + theme_classic()+
		labs(fill = "Reaction ID             ", x ="modelOri", y = "Flux (mmol/gDw/h)")+
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		scale_fill_brewer(palette = "Blues") +
		geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","red"), size=1)

# PIPLC
	ss.data = Test4HumanAstro[Test4HumanAstro$TestSolutionName == 'PIPLC',]						
	mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
            measure.vars=c("Primary", "iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
	mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))     
    mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl',]
	FluxCtrl = mdat_ctrl$value               

	ggplot(mdat, aes(x = variable, y = value)) +
		geom_col(aes(fill = TestSolutionName), width = 0.7, colour = "black") + theme_classic()+
		labs(fill = "Reaction ID             ", x ="modelOri", y = "Flux (mmol/gDw/h)")+
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		scale_fill_brewer(palette = "Blues") +
		geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","red"), size=1)

# PGMT
	ss.data = Test4HumanAstro[Test4HumanAstro$TestSolutionName == 'PGMT',]						
	mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
            measure.vars=c("Primary", "iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
	mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))     
    mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl',]
	FluxCtrl = mdat_ctrl$value               

	ggplot(mdat, aes(x = variable, y = value)) +
		geom_col(aes(fill = TestSolutionName), width = 0.7, colour = "black") + theme_classic()+
		labs(fill = "Reaction ID             ", x ="modelOri", y = "Flux (mmol/gDw/h)")+
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		scale_fill_brewer(palette = "Blues") +
		geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","red"), size=1)

# HMGCOASim
	ss.data = Test4HumanAstro[Test4HumanAstro$TestSolutionName == 'HMGCOASim',]						
	mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
            measure.vars=c("Primary", "iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
	mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))     
    mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl',]
	FluxCtrl = mdat_ctrl$value               

	ggplot(mdat, aes(x = variable, y = value)) +
		geom_col(aes(fill = TestSolutionName), width = 0.7, colour = "black") + theme_classic()+
		labs(fill = "Reaction ID             ", x ="modelOri", y = "Flux (mmol/gDw/h)")+
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		scale_fill_brewer(palette = "Blues") +
		geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","red"), size=1)
