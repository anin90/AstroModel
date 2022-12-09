#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
#~ install.packages("gridExtra")
#~ install.packages("patchwork")
pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table, gridExtra,
patchwork)

setwd("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/5.generateFigures")

pdf("1.ManuscriptFigs/generateFigures_ms.pdf")

#################################
# Load data
#################################
#disrupted in >0 modules/phenotype-of-interest
# Fig.1a
		mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd.csv", header = T, sep = "\t")
		colnames(mat) = c("subSystem", "FVA_BD_abs", "FVA_BD_norm_t1", "FVA_BD_norm_t2", 
							"MTA_BD_abs", "MTA_BD_norm_t1", "MTA_BD_norm_t2")
		mat = subset(mat, TRUE, c("subSystem", "FVA_BD_abs", "FVA_BD_norm_t1", "FVA_BD_norm_t2", 
							"MTA_BD_abs", "MTA_BD_norm_t1", "MTA_BD_norm_t2"))
		row.names(mat) <- mat$subSystem		
		mat = mat[,colnames(mat)!="subSystem"]
		mat = cbind(mat, Total = rowSums(mat!= 0))
		mat = mat[order(mat$Total),]
		keep = rownames(mat)[rowSums(mat)>0]; # filter subSystems disrupted by 0 or more methods
		mat_keep = mat[(row.names(mat) %in% keep),]
		write.table(mat_keep, "1.ManuscriptFigs/data/Tbl_bd_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_r.csv", header = T, sep = "\t")
		colnames(mat) = c("subSystem", "FVA_BD_R_abs", "FVA_BD_R_norm_t1", "FVA_BD_R_norm_t2", 
							"MTA_BD_R_abs", "MTA_BD_R_norm_t1", "MTA_BD_R_norm_t2")
		mat = subset(mat, TRUE, c("subSystem", "FVA_BD_R_abs", "FVA_BD_R_norm_t1", "FVA_BD_R_norm_t2", 
							"MTA_BD_R_abs", "MTA_BD_R_norm_t1", "MTA_BD_R_norm_t2"))
		row.names(mat) <- mat$subSystem
		mat = mat[,colnames(mat)!="subSystem"]
		mat = cbind(mat, Total = rowSums(mat!= 0))
		mat = mat[order(mat$Total),]
		keep = rownames(mat)[rowSums(mat)>0]; # filter subSystems disrupted by 0 or more methods
		mat_keep = mat[(row.names(mat) %in% keep),]
		write.table(mat_keep, "1.ManuscriptFigs/data/Tbl_bd_r_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)	

		mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_nr.csv", header = T, sep = "\t")
		colnames(mat) = c("subSystem", "FVA_BD_NR_abs", "FVA_BD_NR_norm_t1", "FVA_BD_NR_norm_t2", 
							"MTA_BD_NR_abs", "MTA_BD_NR_norm_t1", "MTA_BD_NR_norm_t2")
		mat = subset(mat, TRUE, c("subSystem", "FVA_BD_NR_abs", "FVA_BD_NR_norm_t1", "FVA_BD_NR_norm_t2", 
							"MTA_BD_NR_abs", "MTA_BD_NR_norm_t1", "MTA_BD_NR_norm_t2"))
		row.names(mat) <- mat$subSystem
		mat = mat[,colnames(mat)!="subSystem"]
		mat = cbind(mat, Total = rowSums(mat!= 0))
		mat = mat[order(mat$Total),]
		keep = rownames(mat)[rowSums(mat)>0]; # filter subSystems disrupted by 0 or more methods
		mat_keep = mat[(row.names(mat) %in% keep),]
		write.table(mat_keep, "1.ManuscriptFigs/data/Tbl_bd_nr_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

		BD_tbl <- read.csv("1.ManuscriptFigs/data/Tbl_bd_filt.csv", header = T, sep = "\t")
		BD_R_tbl <- read.csv("1.ManuscriptFigs/data/Tbl_bd_r_filt.csv", header = T, sep = "\t")
		BD_NR_tbl <- read.csv("1.ManuscriptFigs/data/Tbl_bd_nr_filt.csv", header = T, sep = "\t")
		
				# bd	
					df = data.frame("subSystem" = BD_tbl$X, "disruptedBy" = c(BD_tbl$Total))
					bd_lst = rep(df$subSystem, df$disruptedBy)
				# bd_r
					df = data.frame("subSystem" = BD_R_tbl$X, "disruptedBy" = c(BD_R_tbl$Total))
					bd_r_lst = rep(df$subSystem, df$disruptedBy)
				# bd_nr
					df = data.frame("subSystem" = BD_NR_tbl$X, "disruptedBy" = c(BD_NR_tbl$Total))
					bd_nr_lst = rep(df$subSystem, df$disruptedBy)
					
		mat = lst(bd_lst,
					bd_r_lst,
					bd_nr_lst) %>% 	
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)
		write.table(mat, "1.ManuscriptFigs/data/Tbl_mat_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		mat <- read.csv("1.ManuscriptFigs/data/Tbl_mat_filt.csv", header = T, sep = "\t")	
		mat <- mat[, c("name", "bd_lst", "bd_r_lst", "bd_nr_lst")]
		colnames(mat) = c("subSystem", "BD", "BD_R", "BD_NR")
		mm <- melt(mat, id="subSystem")
		plot = ggplot(mm, aes(x=reorder(subSystem, -value), y=value, fill=variable)) + 
			geom_bar(stat="identity", color="black", width = 0.7) + theme_classic() +
			scale_fill_manual(values=c("#339cff", "#33ff58", "#ff4233")) + 
			coord_flip() +
			scale_y_continuous(breaks=c(0:1:10)) + theme(aspect.ratio=1, text = element_text(size=8)) +
			xlab("") + ylab("# of Disrupted Modules") + labs(fill = "Phenotype") +
			theme(axis.text.y=element_text(size=rel(1.1)), legend.position = c(0.7, 0.7)) +
			geom_hline(yintercept=c(2, 3), linetype="dashed", color = c("black","red"), size=0.7)	
			ggsave("1.ManuscriptFigs/Figure1a.pdf", plot, height = 5, width = 7.5)

		plot = ggplot(mm, aes(x=(subSystem), y=value)) + 
			geom_bar(stat="identity", color="black", width = 0.7) + theme_classic()+
			scale_fill_manual(values=c("#339cff", "#33ff58", "#ff4233")) + 
			coord_flip()+
			scale_y_continuous(breaks=c(0:1:10)) +
			xlab("") + ylab("# of Disrupted Modules") + 
			geom_hline(yintercept =2, linetype="dashed", color = c("red"))
		plot + facet_grid(. ~ variable)

#disrupted in >2 modules/phenotype-of-interest
# Fig.1e	
		mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_mat_filt.csv", header = T, sep = "\t")	
		mat <- mat[, c("name", "bd_lst", "bd_r_lst", "bd_nr_lst")]
		colnames(mat) = c("subSystem", "BD", "BD_R", "BD_NR")
		mm2 <- melt(mat, id="subSystem")
		plot = ggplot(mm2, aes(x=reorder(subSystem, -value), y=value, fill=variable)) + 
			geom_bar(stat="identity", color="black", width = 0.7) + theme_classic() +
			scale_fill_manual(values=c("#339cff", "#33ff58", "#ff4233")) + 
			coord_flip() +
			scale_y_continuous(breaks=c(0:1:10)) + theme(aspect.ratio=1) +
			xlab("") + ylab("# of Disrupted Modules") + labs(fill = "Phenotype") +
			theme(axis.text.y=element_text(size=rel(1.1)), legend.position = c(0.8, 0.7)) +
			geom_hline(yintercept=c(2, 3), linetype="dashed", color = c("black","red"), size=0.7)
			ggsave("1.ManuscriptFigs/Figure1e.pdf", plot, height = 2.5, width = 5)
# Fig.1e
		plot = ggplot(mm2, aes(x=(subSystem), y=value)) + 
			geom_bar(stat="identity", color="black", width = 0.7) + theme_classic()+
			scale_fill_manual(values=c("#339cff", "#33ff58", "#ff4233")) + 
			coord_flip()+
			scale_y_continuous(breaks=c(0:1:10)) + theme(aspect.ratio=1) +
			xlab("") + ylab("# of Disrupted Modules") + 
			geom_hline(yintercept =2, linetype="dashed", color = c("red"))
		plot + facet_grid(. ~ variable)

		mm = subset(mm, subSystem %in% mm2$subSystem)
		plot = ggplot(mm, aes(x=(subSystem), y=value)) + 
			geom_bar(stat="identity", color="black", width = 0.7) + theme_classic()+
			scale_fill_manual(values=c("#339cff", "#33ff58", "#ff4233")) + 
			coord_flip()+
			scale_y_continuous(breaks=c(0:1:10)) + theme(aspect.ratio=1) +
			xlab("") + ylab("# of Disrupted Modules") + 
			geom_hline(yintercept =2, linetype="dashed", color = c("red"))
		plot + facet_grid(. ~ variable)
