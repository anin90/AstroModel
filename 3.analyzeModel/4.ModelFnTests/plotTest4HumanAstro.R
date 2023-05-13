#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
#~ install.packages("gridExtra")
#~ install.packages("patchwork")
#~ install.packages("cowplot")
pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table, gridExtra,
patchwork, cowplot)

	# Add a common legend for multiple ggplot2 graphs
			library(gridExtra)
			get_legend<-function(myggplot){
			  tmp <- ggplot_gtable(ggplot_build(myggplot))
			  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
			  legend <- tmp$grobs[[leg]]
			  return(legend)
			}

setwd("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/4.ModelFnTests")

pdf("plotTest4HumanAstro.pdf")

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

	# BD_PI5P4K
		df = Test4HumanAstro_abs
		ss.data = df[df$TestSolutionName == 'PI5P4K',]
		mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
				measure.vars=c("Primary_Ctrl", "iPS_Ctrl_a", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
				measure.vars=c("iPS_Ctrl_a", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl_a',]
		FluxCtrl = mdat_ctrl$value

		BD_PI5P4K= ggplot(mdat, aes(x = variable, y = value)) +
			geom_col(aes(fill = TestSolutionName), width = 0.7, color = 'red') + theme_classic()+
			labs(fill = "Reaction ID", x ="", y = "") +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			scale_fill_brewer(palette = "Blues") +
			# scale_fill_manual(values = c("Primary_Ctrl" = "blue", "iPS_Ctrl" = "blue", "iPS_BD" = "blue", "iPS_BD_R" = "blue",  "iPS_BD_NR" = "red")) +
			geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","black"), size=1) + 
			theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
			theme(legend.position = "none")

	# BD_HMGCOASim
		df = Test4HumanAstro_abs
		ss.data = df[df$TestSolutionName == 'HMGCOASim',]
		mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
				measure.vars=c("Primary_Ctrl", "iPS_Ctrl_a", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
				measure.vars=c("iPS_Ctrl_a", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl_a',]
		FluxCtrl = mdat_ctrl$value

		BD_HMGCOASim = ggplot(mdat, aes(x = variable, y = value)) +
			geom_col(aes(fill = TestSolutionName), width = 0.7, color = 'red') + theme_classic()+
			labs(fill = "Reaction ID", x ="", y = "") +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			scale_fill_brewer(palette = "Blues") +
			# scale_fill_manual(values = c("Primary_Ctrl" = "blue", "iPS_Ctrl" = "blue", "iPS_BD" = "blue", "iPS_BD_R" = "blue",  "iPS_BD_NR" = "red")) +
			geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","black"), size=1) + 
			theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
			theme(legend.position = "none")

	# BD_FAOXC164C165m
		df = Test4HumanAstro_norm_t2
		ss.data = df[df$TestSolutionName == 'FAOXC164C165m',]
		mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
				measure.vars=c("Primary_Ctrl", "iPS_Ctrl_a", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
				measure.vars=c("iPS_Ctrl_a", "iPS_BD", "iPS_BD_R", "iPS_BD_NR"))
		mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl_a',]
		FluxCtrl = mdat_ctrl$value

		BD_FAOXC164C165m = ggplot(mdat, aes(x = variable, y = value)) +
			geom_col(aes(fill = TestSolutionName), width = 0.7, color = 'red') + theme_classic()+
			labs(fill = "Reaction ID", x ="", y = "") +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			scale_fill_brewer(palette = "Blues") +
			# scale_fill_manual(values = c("Primary_Ctrl" = "blue", "iPS_Ctrl" = "blue", "iPS_BD" = "blue", "iPS_BD_R" = "blue",  "iPS_BD_NR" = "red")) +
			geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","black"), size=1) + 
			theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
			theme(legend.position = "none")

	# ST_FAOXC164C165m
		df = Test4HumanAstro_norm_t2
		ss.data = df[df$TestSolutionName == 'FAOXC164C165m',]
		mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
				measure.vars=c("Primary_Ctrl", "iPS_Ctrl_b", "iPS_HT", "iPS_ST"))
		mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
				measure.vars=c("iPS_Ctrl_b", "iPS_HT", "iPS_ST"))
		mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl_b',]
		FluxCtrl = mdat_ctrl$value

		ST_FAOXC164C165m = ggplot(mdat, aes(x = variable, y = value)) +
			geom_col(aes(fill = TestSolutionName), width = 0.7, color = 'red') + theme_classic()+
			labs(fill = "Reaction ID", x ="", y = "") +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			scale_fill_brewer(palette = "Blues") +
			# scale_fill_manual(values = c("Primary_Ctrl" = "blue", "iPS_Ctrl" = "blue", "iPS_BD" = "blue", "iPS_BD_R" = "blue",  "iPS_BD_NR" = "red")) +
			geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","black"), size=1) + 
			theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
			theme(legend.position = "none")

	# ST_LNLCCPT2
		df = Test4HumanAstro_norm_t1
		ss.data = df[df$TestSolutionName == 'LNLCCPT2',]
		mdat = melt(ss.data, id.vars=c("TestSolutionName", "TestSolutionGroup"),
				measure.vars=c("Primary_Ctrl", "iPS_Ctrl_b", "iPS_HT", "iPS_ST"))
		mdat_ctrl = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
				measure.vars=c("iPS_Ctrl_b", "iPS_HT", "iPS_ST"))
		mdat_ctrl = mdat_ctrl[mdat_ctrl$variable == 'iPS_Ctrl_b',]
		FluxCtrl = mdat_ctrl$value

		ST_LNLCCPT2 = ggplot(mdat, aes(x = variable, y = value)) +
			geom_col(aes(fill = TestSolutionName), width = 0.7, color = 'red') + theme_classic()+
			labs(fill = "Reaction ID", x ="", y = "") +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			scale_fill_brewer(palette = "Blues") +
			# scale_fill_manual(values = c("Primary_Ctrl" = "blue", "iPS_Ctrl" = "blue", "iPS_BD" = "blue", "iPS_BD_R" = "blue",  "iPS_BD_NR" = "red")) +
			geom_hline(yintercept=c(FluxCtrl/1.5, FluxCtrl/0.8), linetype="dashed", color = c("black","black"), size=1) + 
			theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
			theme(legend.position = "none")

#############################################
# print plots - grid/ stitch/ patchwork/ etc..
#############################################

				# plotGrid
						
						plot_grid(BD_PI5P4K, BD_HMGCOASim, BD_FAOXC164C165m, 
										ST_FAOXC164C165m, ST_LNLCCPT2,
										labels=c("a.) BD_PI5P4K", "b.) BD_HMGCOASim", "c.) BD_FAOXC164C165m", 
														"d.) SCZ_FAOXC164C165m", "e.) SCZ_LNLCCPT2" ), 
										label_size = 6, ncol = 4, nrow = 2, label_x = 0, label_y = 0, hjust = -0.5, vjust = -0.5)
										
						svg("plotTest4HumanAstroFigs/1a.svg")
						plot_grid(BD_PI5P4K, BD_HMGCOASim, BD_FAOXC164C165m, 
										ST_FAOXC164C165m, ST_LNLCCPT2,
										labels=c("a.) BD_PI5P4K", "b.) BD_HMGCOASim", "c.) BD_FAOXC164C165m", 
														"d.) SCZ_FAOXC164C165m", "e.) SCZ_LNLCCPT2" ), 
										label_size = 6, ncol = 4, nrow = 2, label_x = 0, label_y = 0, hjust = -0.5, vjust = -0.5)
						dev.off()
