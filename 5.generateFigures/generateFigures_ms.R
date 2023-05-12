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

setwd("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/5.generateFigures")

pdf("1.ManuscriptFigs/generateFigures_ms.pdf")

#################################
# Load data
#################################
#disrupted in >0 modules/phenotype-of-interest
# Fig.1a
		mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd.csv", header = T, sep = "\t")
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

		mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_r.csv", header = T, sep = "\t")
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

		mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_nr.csv", header = T, sep = "\t")
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

		svg("1.ManuscriptFigs/1a.svg")
		plot + facet_grid(. ~ variable)
		dev.off() 

#disrupted in >2 modules/phenotype-of-interest
# Fig.1e	
		mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_mat_filt.csv", header = T, sep = "\t")	
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
		
		svg("1.ManuscriptFigs/1e.svg")
		plot + facet_grid(. ~ variable)
		dev.off() 

#Rxns belonging to the subSystems-disrupted in >2 modules/phenotype-of-interest
		BD <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/BD_Rxns.csv", header = T, sep = "\t")
		BD_R <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/BD_R_Rxns.csv", header = T, sep = "\t")
		BD_NR <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/BD_NR_Rxns.csv", header = T, sep = "\t")

		head(BD)		
		dim(BD)
		dim(BD_R)
		dim(BD_NR)

###########################
# Fig.2 - compare model stats
###########################
# Load data
		modelStats_pre = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/modelStatsMatFiles/modelStatsMatSol.csv", header = T, sep = "\t")
		modelStats_post = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/4.modelComparison/modelStatsMatFilesFinal/modelStatsMatSolFinal.csv", header = T, sep = "\t")
		modelStats_predecessor = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/4.modelComparison/predecessorModel/predecessorModelComparison.csv", header = T, sep = "\t")

# Fig.2 - test
		ggplot(modelStats_post, aes(x=Phenotype, y=fluxInconsistentRxnsPrct, fill=Dataset)) + 
			geom_bar(stat="identity", position=position_dodge(), width = 0.5, color="black")+
			scale_fill_manual(values=c('#999999','#E69F00', '#56B4E9')) +
			theme_classic() + theme(aspect.ratio=1) + xlab("") + ylab("fluxInconsistentRxnsPrct") +
			coord_flip() + facet_grid(. ~ ExpThreshold)
			
		ggplot(modelStats_post, aes(x=Phenotype, y=fluxInconsistentRxnsPrct, fill=ExpThreshold)) + 
			geom_bar(stat="identity", position=position_dodge(), width = 0.5, color="black")+
			scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
			theme_classic() + theme(aspect.ratio=1) + xlab("") + ylab("fluxInconsistentRxnsPrct") +
			coord_flip()

# Fig.2a-e ------Compare models pre & post expansion------

		# merge colNames to create modelID (pre expansion)
		pre = modelStats_pre
		pre$Model_ID <- paste(pre$Dataset, pre$Phenotype, pre$ExpThreshold, pre$MEM, sep="_")
		pre <- pre %>% relocate(Model_ID, .before = Dataset)
		
		# merge colNames to create modelID (post expansion)
		post = modelStats_post
		post$Model_ID <- paste(post$Dataset, post$Phenotype, post$ExpThreshold, post$MEM, sep="_")
		post <- post %>% relocate(Model_ID, .before = Dataset)

		# merge pre- & post- exansion modelStats
		dat = merge(pre,post, by.x = "Model_ID", by.y = "Model_ID")

				# plot the "Number of modelRxns" - pre & post
				mat <- dat[, c("Model_ID", "Dataset.x", "Phenotype.x", "ExpThreshold.x", "MEM.x", "modelRxns.x", "modelRxns.y")]
				colnames(mat) = c("Model_ID", "Dataset", "Phenotype", "ExpThreshold", "MEM", "preExpansion", "postExpansion")
				mat = mat %>% pivot_longer(cols=c('preExpansion', 'postExpansion'), names_to='stage', values_to='modelRxns')
				
						a1 = ggplot(mat, aes(x=stage, y=modelRxns, fill = Phenotype)) +
							geom_boxplot(position=position_dodge(0.5), width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("# of rxns in model") + 
							theme(text = element_text(size = 18)) +					
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()
						
						a2 = ggplot(mat, aes(x=stage, y=modelRxns, fill = Phenotype)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("# of rxns in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()
							
						a = ggplot(mat, aes(x=stage, y=modelRxns, fill = Phenotype)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_point(aes(shape=ExpThreshold, color=Phenotype), size=4, position=position_dodge(0.5)) +
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("# of rxns in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/2a.svg")
						a
						dev.off()

				# plot the "Number of modelMets" - pre & post
				mat <- dat[, c("Model_ID", "Dataset.x", "Phenotype.x", "ExpThreshold.x", "MEM.x", "modelMets.x", "modelMets.y")]
				colnames(mat) = c("Model_ID", "Dataset", "Phenotype", "ExpThreshold", "MEM", "preExpansion", "postExpansion")
				mat = mat %>% pivot_longer(cols=c('preExpansion', 'postExpansion'), names_to='stage', values_to='modelMets')

						b = ggplot(mat, aes(x=stage, y=modelMets, fill = Phenotype)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_point(aes(shape=ExpThreshold, color=Phenotype), size=4, position=position_dodge(0.5)) +
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("# of mets in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							scale_y_continuous(breaks = seq(3800, 4200, by = 200)) + 
							coord_flip()

						svg("1.ManuscriptFigs/2b.svg")
						b
						dev.off()

				# plot the "Number of fluxInconsistentRxnsPrct" - pre & post
				mat <- dat[, c("Model_ID", "Dataset.x", "Phenotype.x", "ExpThreshold.x", "MEM.x", "fluxInconsistentRxnsPrct.x", "fluxInconsistentRxnsPrct.y")]
				colnames(mat) = c("Model_ID", "Dataset", "Phenotype", "ExpThreshold", "MEM", "preExpansion", "postExpansion")
				mat = mat %>% pivot_longer(cols=c('preExpansion', 'postExpansion'), names_to='stage', values_to='fluxInconsistentRxnsPrct')

						c = ggplot(mat, aes(x=stage, y=fluxInconsistentRxnsPrct, fill = Phenotype)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_point(aes(shape=ExpThreshold, color=Phenotype), size=4, position=position_dodge(0.5)) +
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("% of fluxInconsistent rxns in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/2c.svg")
						c
						dev.off()

				# plot the "Number of overlapCoreRxnsPrct" - pre & post
				mat <- dat[, c("Model_ID", "Dataset.x", "Phenotype.x", "ExpThreshold.x", "MEM.x", "overlapCoreRxnsPrct.x", "overlapCoreRxnsPrct.y")]
				colnames(mat) = c("Model_ID", "Dataset", "Phenotype", "ExpThreshold", "MEM", "preExpansion", "postExpansion")
				mat = mat %>% pivot_longer(cols=c('preExpansion', 'postExpansion'), names_to='stage', values_to='overlapCoreRxnsPrct')

						d = ggplot(mat, aes(x=stage, y=overlapCoreRxnsPrct, fill = Phenotype)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_point(aes(shape=ExpThreshold, color=Phenotype), size=4, position=position_dodge(0.5)) +
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("% of core rxns in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/2d.svg")
						d
						dev.off()

				# plot the "Number of overlapLewisPrct" - pre & post
				mat <- dat[, c("Model_ID", "Dataset.x", "Phenotype.x", "ExpThreshold.x", "MEM.x", "overlapLewisPrct.x", "overlapLewisPrct.y")]
				colnames(mat) = c("Model_ID", "Dataset", "Phenotype", "ExpThreshold", "MEM", "preExpansion", "postExpansion")
				mat = mat %>% pivot_longer(cols=c('preExpansion', 'postExpansion'), names_to='stage', values_to='overlapLewisPrct')

						e = ggplot(mat, aes(x=stage, y=overlapLewisPrct, fill = Phenotype)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_point(aes(shape=ExpThreshold, color=Phenotype), size=4, position=position_dodge(0.5)) +
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("% of astrocytic rxns, from Lewis et al. 2010, in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/2e.svg")
						e
						dev.off()

# Fig.2f-g ------Compare models (postExpansion) with predecessors------

				# plot the "Number of modelRxns" - only post
				mat <- dat[, c("Model_ID", "Dataset.x", "Phenotype.x", "ExpThreshold.x", "MEM.x", "modelRxns.y")]
				colnames(mat) = c("Model_ID", "Dataset", "Phenotype", "ExpThreshold", "MEM", "postExpansion")
				mat_Rxns = mat %>% pivot_longer(cols=c('postExpansion'), names_to='stage', values_to='modelRxns')
				
				# plot the "Number of modelMets" - only post
				mat <- dat[, c("Model_ID", "Dataset.x", "Phenotype.x", "ExpThreshold.x", "MEM.x", "modelMets.y")]
				colnames(mat) = c("Model_ID", "Dataset", "Phenotype", "ExpThreshold", "MEM", "postExpansion")
				mat_Mets = mat %>% pivot_longer(cols=c('postExpansion'), names_to='stage', values_to='modelMets')
				
				# compare modelStats of "Chellappa" vs "Predecessors"
				mat = modelStats_predecessor %>% 
																		add_row(Model_Source = rep(c('Chellappa. (2023)'), nrow(mat_Rxns)),
																		Model_Source_Anno = rep(c('This Study'), nrow(mat_Rxns)),
																		Phenotype = mat_Rxns$"Phenotype",
																		modelRxns = mat_Rxns$"modelRxns",
																		modelMets = mat_Mets$"modelMets")

						# of rxns														
						f = ggplot(mat, aes(x=Model_Source_Anno, y=modelRxns, fill = Model_Source)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							scale_fill_discrete(breaks=c('Chellappa. (2023)', 'Baloni. (2020)', 'Martín-Jiménez. (2017)', 'DiNuzzo. (2017)', 'Lewis. (2010)', 'Çakir. (2007)')) +
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("# of rxns in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							scale_y_continuous(breaks = seq(0, 6000, by = 2000)) + 
							coord_flip()

						svg("1.ManuscriptFigs/2f.svg")
						f
						dev.off()
																				
						# of mets														
						g = ggplot(mat, aes(x=Model_Source_Anno, y=modelMets, fill = Model_Source)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							scale_fill_discrete(breaks=c('Chellappa. (2023)', 'Baloni. (2020)', 'Martín-Jiménez. (2017)', 'DiNuzzo. (2017)', 'Lewis. (2010)', 'Çakir. (2007)')) +
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("# of mets in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							scale_y_continuous(breaks = seq(0, 6000, by = 2000)) + 
							coord_flip()
							
						svg("1.ManuscriptFigs/2g.svg")
						g
						dev.off()

##########################################
# Fig.3 ------Compare models - Koskuvi. (2022) vs Biju. (2020)------
##########################################

					# merge colNames to create modelID (pre expansion)
					pre = modelStats_pre
					pre$Model_ID <- paste(pre$Dataset, pre$Phenotype, pre$ExpThreshold, pre$MEM, sep="_")
					pre <- pre %>% relocate(Model_ID, .before = Dataset)

					# subset only Koskuvi and Biju models (pre expansion)
					pre = pre %>% filter(grepl('Koskuvi|Biju', Model_ID))
					mat = pre
					
						# of rxns
						h = ggplot(mat, aes(x=Dataset, y=modelRxns, fill = ExpThreshold)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("# of rxns in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/3a.svg")
						h
						dev.off()

						# of mets
						i = ggplot(mat, aes(x=Dataset, y=modelMets, fill = ExpThreshold)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("# of mets in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/3b.svg")
						i
						dev.off()

						# % of fluxInconsistent rxns in model
						j = ggplot(mat, aes(x=Dataset, y=fluxInconsistentRxnsPrct, fill = ExpThreshold)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("% of fluxInconsistent rxns in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/3c.svg")
						j
						dev.off()

						# % of fluxInconsistent rxns in model
						k = ggplot(mat, aes(x=Dataset, y=overlapCoreRxnsPrct, fill = ExpThreshold)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("% of core rxns in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/3d.svg")
						k
						dev.off()

						# % of fluxInconsistent rxns in model
						l = ggplot(mat, aes(x=Dataset, y=overlapLewisPrct, fill = ExpThreshold)) +
							geom_boxplot(fill="white", width = 0.5, color="black")+
							geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.5)) + 
							theme_classic() + theme(aspect.ratio=0.5) + xlab("") + ylab("% of astrocytic rxns, from Lewis et al. 2010, in model") + 
							theme(text = element_text(size = 18)) +
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							coord_flip()

						svg("1.ManuscriptFigs/3e.svg")
						l
						dev.off()

#############################################
# Fig.4 ------Merge final results - Vadodaria (BD) and Koskuvi. (SCZ)------
#############################################

		# load table
		BD <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_subSystem_pval.csv", header = T, sep = "\t")
		BD_R <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_subSystem_pval.csv", header = T, sep = "\t")
		BD_NR <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_subSystem_pval.csv", header = T, sep = "\t")
		ST <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/PlotResults/plotDataOverlaysHyper_Tbl_Final/ST_subSystem_pval.csv", header = T, sep = "\t")

		dat = rbind(BD, BD_R, BD_NR, ST)
		dat$SignificantBy = as.character(dat$No.of.disrupted.modules)
		dat$SignificantBy <- paste0(dat$SignificantBy, " modules")
		dat

						# grid plot - merge vadodaria & koskuvi results
						m = ggplot(dat, aes(x=(subSystem), y=-log10(rowMeans),  fill = SignificantBy)) + 
							geom_bar(stat="identity", color="black", width = 0.7) + theme_classic()+
							scale_fill_manual(values=c('#999999','#E69F00', '#56B4E9')) + 
							coord_flip()+
							xlab("") + ylab("-log10(mean.fdr.p.val)") + 
							geom_hline(yintercept = -log10(0.05), linetype="dashed", color = c("red")) + 
							theme(text=element_text(color="black"), axis.text=element_text(color="black")) +
							facet_grid(. ~ Phenotype)
						
						svg("1.ManuscriptFigs/4a.svg")
						m + theme(aspect.ratio=2, text = element_text(size=7)) +
								theme(axis.text.y=element_text(size=rel(1.4)), axis.text.x=element_text(size=rel(1.4))) +
								theme(legend.text=element_text(size=7)) +
								theme(strip.text.x = element_text(size = 7)) +
								theme(axis.title.y=element_text(size=7))
						dev.off()

#############################################
# print plots - grid/ stitch/ patchwork/ etc..
#############################################

				# plotGrid
				plot_grid(a1, a2, labels=c("a1.)", "a2.)"), ncol = 1, nrow = 2)
				plot_grid(a, b, labels=c("2a.)", "2b.)"), ncol = 1, nrow = 2)
				plot_grid(c, d, labels=c("2c.)", "2d.)"), ncol = 1, nrow = 2)
				plot_grid(e, labels=c("2e.)"), ncol = 1, nrow = 2)
				plot_grid(f, g, labels=c("2f.)", "2g.)"), ncol = 1, nrow = 2)
				plot_grid(h, i, labels=c("3a.)", "3b.)"), ncol = 1, nrow = 2)
				plot_grid(j, k, labels=c("3c.)", "3d.)"), ncol = 1, nrow = 2)
				plot_grid(l, labels=c("3e.)"), ncol = 1, nrow = 2)
				plot_grid(m, labels=c("4a.)"), ncol = 1, nrow = 2)
				
				# patchwork
				x1 =	 (a +theme(legend.position="none")) / 
							(b+theme(legend.position="none")) | 
							((c +theme(legend.position="none")) / 
							(d+theme(legend.position="none"))) + plot_layout(guides = "collect")

				svg("1.ManuscriptFigs/Fig2.svg")
				x1
				dev.off()
