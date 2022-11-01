#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table) 

setwd("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/4.ModelFnTests")

pdf("plotTest4HumanAstro.pdf", width=3.5, height=3.5)

######################
# Load data - TestResults_files
######################

	# load data
		Test4HumanAstro_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/4.ModelFnTests/TestResults_files/Test4HumanAstro_abs.csv", header = T, sep = "\t")
		Test4HumanAstro_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/4.ModelFnTests/TestResults_files/Test4HumanAstro_norm_t1.csv", header = T, sep = "\t")
		Test4HumanAstro_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/4.ModelFnTests/TestResults_files/Test4HumanAstro_norm_t2.csv", header = T, sep = "\t")

######################
# Plot data - TestResults_files
######################

	# PI5P4K
		df = Test4HumanAstro_abs
		ss.data = df[df$TestSolutionName == 'PI5P4K',]
		mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
				measure.vars=c("Primary_Ctrl", "iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
				measure.vars=c("iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl',]
		FluxCtrl = mdat_ctrl$value

		ggplot(mdat, aes(x = variable, y = value)) +
			geom_col(aes(fill = TestSolutionName), width = 0.7, color = 'red') + theme_classic()+
			labs(fill = "Reaction ID", x ="", y = "") +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			scale_fill_brewer(palette = "Blues") +
			# scale_fill_manual(values = c("Primary_Ctrl" = "blue", "iPS_Ctrl" = "blue", "iPS_BD" = "blue", "iPS_BD_R" = "blue",  "iPS_BD_NR" = "red")) +
			geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","black"), size=1)

	# HMGCOASim
		df = Test4HumanAstro_abs
		ss.data = df[df$TestSolutionName == 'HMGCOASim',]
		mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
				measure.vars=c("Primary_Ctrl", "iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
				measure.vars=c("iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl',]
		FluxCtrl = mdat_ctrl$value

		ggplot(mdat, aes(x = variable, y = value)) +
			geom_col(aes(fill = TestSolutionName), width = 0.7, color = 'red') + theme_classic()+
			labs(fill = "Reaction ID", x ="", y = "") +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			scale_fill_brewer(palette = "Blues") +
			# scale_fill_manual(values = c("Primary_Ctrl" = "blue", "iPS_Ctrl" = "blue", "iPS_BD" = "blue", "iPS_BD_R" = "blue",  "iPS_BD_NR" = "red")) +
			geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","black"), size=1)

	# FAOXC164C165m
		df = Test4HumanAstro_norm_t2
		ss.data = df[df$TestSolutionName == 'FAOXC164C165m',]
		mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
				measure.vars=c("Primary_Ctrl", "iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
				measure.vars=c("iPS_Ctrl", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl',]
		FluxCtrl = mdat_ctrl$value

		ggplot(mdat, aes(x = variable, y = value)) +
			geom_col(aes(fill = TestSolutionName), width = 0.7, color = 'red') + theme_classic()+
			labs(fill = "Reaction ID", x ="", y = "") +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			scale_fill_brewer(palette = "Blues") +
			# scale_fill_manual(values = c("Primary_Ctrl" = "blue", "iPS_Ctrl" = "blue", "iPS_BD" = "blue", "iPS_BD_R" = "blue",  "iPS_BD_NR" = "red")) +
			geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","black"), size=1)





