#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table) 

setwd("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria")

pdf("PlotResults/filterDataOverlays.pdf")

#################################
# Load tbl - plotDataOverlaysHyper_Tbl_Final
#################################

	# load table
		BD_S1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S1.csv", header = T, sep = "\t")
		BD_S2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S2.csv", header = T, sep = "\t")
		BD_S3 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S3.csv", header = T, sep = "\t")
		BD_R_S1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S1.csv", header = T, sep = "\t")
		BD_R_S2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S2.csv", header = T, sep = "\t")
		BD_R_S3 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S3.csv", header = T, sep = "\t")
		BD_NR_S1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S1.csv", header = T, sep = "\t")
		BD_NR_S2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S2.csv", header = T, sep = "\t")
		BD_NR_S3 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S3.csv", header = T, sep = "\t")
		BD_NR_S4 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S4.csv", header = T, sep = "\t")
		BD_NR_S5 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S5.csv", header = T, sep = "\t")
		BD_NR_S6 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S6.csv", header = T, sep = "\t")
		BD_NR_S7 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S7.csv", header = T, sep = "\t")
		BD_NR_S8 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S8.csv", header = T, sep = "\t")

#################################
# Plot tbl - plotDataOverlaysHyper_Tbl_Final
#################################

	# BD

		# Methionine and cysteine metabolism (BD_S1)
			df = as.data.frame(table(BD_S1$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Methionine and cysteine metabolism (BD_S1)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_S1_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		# Keratan sulfate degradation (BD_S2)
			df = as.data.frame(table(BD_S2$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Keratan sulfate degradation (BD_S2)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_S2_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		# Fatty acid oxidation (BD_S3)
			df = as.data.frame(table(BD_S3$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Fatty acid oxidation (BD_S3)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_S3_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
				
	# BD_R

		# ROS detoxification (BD_R_S1)
			df = as.data.frame(table(BD_R_S1$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="ROS detoxification (BD_R_S1)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_R_S1_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			
		# Drug metabolism (BD_R_S2)
			df = as.data.frame(table(BD_R_S2$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Drug metabolism (BD_R_S2)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_R_S2_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			
		# Transport, extracellular (BD_R_S3)
			df = as.data.frame(table(BD_R_S3$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Transport, extracellular (BD_R_S3)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_R_S3_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

	# BD_NR

		# Cholesterol metabolism (BD_NR_S1)
			df = as.data.frame(table(BD_NR_S1$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Cholesterol metabolism (BD_NR_S1)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_NR_S1_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			
		# Inositol phosphate metabolism (BD_NR_S2)
			df = as.data.frame(table(BD_NR_S2$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Inositol phosphate metabolism (BD_NR_S2)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_NR_S2_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		# Fatty acid oxidation (BD_NR_S3)
			df = as.data.frame(table(BD_NR_S3$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Fatty acid oxidation (BD_NR_S3)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_NR_S3_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		# Heparan sulfate degradation (BD_NR_S4)
			df = as.data.frame(table(BD_NR_S4$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Heparan sulfate degradation (BD_NR_S4)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_NR_S4_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		# Tyrosine metabolism (BD_NR_S5)
			df = as.data.frame(table(BD_NR_S5$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Tyrosine metabolism (BD_NR_S5)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_NR_S5_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		# N-glycan synthesis (BD_NR_S6)
			df = as.data.frame(table(BD_NR_S6$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="N-glycan synthesis (BD_NR_S6)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_NR_S6_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		# Phosphatidylinositol phosphate metabolism (BD_NR_S7)
			df = as.data.frame(table(BD_NR_S7$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Phosphatidylinositol phosphate metabolism (BD_NR_S7)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_NR_S7_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		# Keratan sulfate degradation (BD_NR_S8)
			df = as.data.frame(table(BD_NR_S8$rxnList))
			ggplot(df, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="identity", aes(fill=Freq>1), width = 0.7) + theme_classic() +
				coord_flip() + scale_y_continuous(breaks=c(0:1:6)) + scale_fill_manual(values = c('#BBBBBC', '#464546')) +
				theme(axis.text.y=element_text(size=rel(1.1))) + labs(title="Keratan sulfate degradation (BD_NR_S8)", x="", y = "# of Disrupted Modules")
			df_desc = df[order(df$Freq, decreasing = TRUE), ]
			write.table(df_desc, "PlotResults/filterDataOverlays_Tbl/BD_NR_S8_desc.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)



