#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table) 

setwd("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi")

pdf("PlotResults/plotDataOverlaysHyper.pdf")

#################################################################################################
# https://stats.stackexchange.com/questions/95523/significance-of-overlap-between-multiple-lists
#################################################################################################
	hyper_matrix <- function(gene.list, background){
		
	  # generate every combinations of two gene lists
	  combination <- expand.grid(names(gene.list),names(gene.list))
	  combination$values <- rep(NA, times=nrow(combination))

	  # convert long table into wide
	  combination <- reshape(combination, idvar="Var1", timevar="Var2", direction="wide")
	  rownames(combination) <- combination$Var1
	  combination <- combination[,-1]
	  colnames(combination) <- gsub("values.", "", colnames(combination))

	  # calculate the length of overlap of each pair
	  for(i in colnames(combination)){
		for(j in rownames(combination)){
		  combination[j,i]<-length(intersect(gene.list[[j]],gene.list[[i]]))
		}
	  }

	  # calculate the significance of the overlap of each pair
	  for(m in 1:length(gene.list)){
		for(n in 1:length(gene.list)){
		  if(n>m){
			combination[n,m] <- phyper(combination[m,n]-1, length(gene.list[[m]]), background-length(gene.list[[m]]), length(gene.list[[n]]), lower.tail=F)
			
			# note that the phyper function (lower.tail=F) give the probability of P[X>x], so the the overlap length should subtract 1 to get a P[X>=x].
		  }
		}
	  }
	  return(combination)
	}


		RBIND <- function(datalist) {
					  require(plyr)
					  temp <- rbind.fill(datalist)
					  temp
					}

###################
# Load rxn data
###################

	# model_abs
		# background
		All_Primary_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_Primary_TP_abs.csv", header = T, sep = "\t")
		All_Ctrl_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_Ctrl_TP_abs.csv", header = T, sep = "\t")
		All_HT_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_HT_TP_abs.csv", header = T, sep = "\t")
		All_ST_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_ST_TP_abs.csv", header = T, sep = "\t")
		# results_koskuvi
		FVA_ST_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_significant_abs/ST_rxns_fdr_abs.csv", header = T, sep = "\t")
		MTA_ST_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_prctile_top_significant_abs/ST_rxns_fdr.csv", header = T, sep = "\t")
		# results_lithium
		Lanz_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Lanz_vs_Ctrl_abs.csv", header = T, sep = "\t")
		Rivera_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Rivera_vs_Ctrl_abs.csv", header = T, sep = "\t")
		Akkouh_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Akkouh_vs_Ctrl_abs.csv", header = T, sep = "\t")
		
	# model_norm_t1
		# background
		All_Primary_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_Primary_TP_norm_t1.csv", header = T, sep = "\t")
		All_Ctrl_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_Ctrl_TP_norm_t1.csv", header = T, sep = "\t")
		All_HT_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_HT_TP_norm_t1.csv", header = T, sep = "\t")
		All_ST_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_ST_TP_norm_t1.csv", header = T, sep = "\t")
		# results_koskuvi
		FVA_ST_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_significant_norm_t1/ST_rxns_fdr_norm_t1.csv", header = T, sep = "\t")
		MTA_ST_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_prctile_top_significant_norm_t1/ST_rxns_fdr.csv", header = T, sep = "\t")
		# results_lithium
		Lanz_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Lanz_vs_Ctrl_norm_t1.csv", header = T, sep = "\t")
		Rivera_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Rivera_vs_Ctrl_norm_t1.csv", header = T, sep = "\t")
		Akkouh_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Akkouh_vs_Ctrl_norm_t1.csv", header = T, sep = "\t")

	# model_norm_t2
		# background
		All_Primary_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_Primary_TP_norm_t2.csv", header = T, sep = "\t")
		All_Ctrl_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_Ctrl_TP_norm_t2.csv", header = T, sep = "\t")
		All_HT_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_HT_TP_norm_t2.csv", header = T, sep = "\t")
		All_ST_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_ST_TP_norm_t2.csv", header = T, sep = "\t")
		# results_koskuvi
		FVA_ST_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_significant_norm_t2/ST_rxns_fdr_norm_t2.csv", header = T, sep = "\t")
		MTA_ST_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_prctile_top_significant_norm_t2/ST_rxns_fdr.csv", header = T, sep = "\t")
		# results_lithium
		Lanz_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Lanz_vs_Ctrl_norm_t2.csv", header = T, sep = "\t")
		Rivera_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Rivera_vs_Ctrl_norm_t2.csv", header = T, sep = "\t")
		Akkouh_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/AnalysisFiles/Akkouh_vs_Ctrl_norm_t2.csv", header = T, sep = "\t")

############################
# ST vs Lithium
############################
	
	# model_abs

		x = list(Lanz = Lanz_abs$Var1, Rivera = Rivera_abs$Var1, Akkouh = Akkouh_abs$Var1,
					FVA_ST = FVA_ST_abs$rxnList, MTA_ST = MTA_ST_abs$del_rxnID_KO)

		y = list(All_Primary_abs$Rxn, All_Ctrl_abs$Rxn, All_HT_abs$Rxn, All_ST_abs$Rxn)
		yInt = Reduce(intersect,y); background = length(yInt);
		
		M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
		M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
		rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
		pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
		main = 'Abs - ST vs Li+', #color=colorRampPalette(c("green", "yellow", "orangered"))(50))
		color = c("green", "orange"), breaks = c(0, 0.05, 1), legend = F)

		M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;
		
	# model_norm_t1

		x = list(Lanz = Lanz_norm_t1$Var1, Rivera = Rivera_norm_t1$Var1, Akkouh = Akkouh_norm_t1$Var1,
					FVA_ST = FVA_ST_norm_t1$rxnList, MTA_ST = MTA_ST_norm_t1$del_rxnID_KO)

		y = list(All_Primary_norm_t1$Rxn, All_Ctrl_norm_t1$Rxn, All_HT_norm_t1$Rxn, All_ST_norm_t1$Rxn)
		yInt = Reduce(intersect,y); background = length(yInt);
		
		M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
		M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
		rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
		pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
		main = 'Norm_t1 - ST vs Li+', #color=colorRampPalette(c("green", "yellow", "orangered"))(50))
		color = c("green", "orange"), breaks = c(0, 0.05, 1), legend = F)

		M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;

	# model_norm_t2

		x = list(Lanz = Lanz_norm_t2$Var1, Rivera = Rivera_norm_t2$Var1, Akkouh = Akkouh_norm_t2$Var1,
					FVA_ST = FVA_ST_norm_t2$rxnList, MTA_ST = MTA_ST_norm_t2$del_rxnID_KO)

		y = list(All_Primary_norm_t2$Rxn, All_Ctrl_norm_t2$Rxn, All_HT_norm_t2$Rxn, All_ST_norm_t2$Rxn)
		yInt = Reduce(intersect,y); background = length(yInt);
		
		M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
		M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
		rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
		pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
		main = 'Norm_t2 - ST vs Li+', #color=colorRampPalette(c("green", "yellow", "orangered"))(50))
		color = c("green", "orange"), breaks = c(0, 0.05, 1), legend = F)

		M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;

	# model_all

		x = list(Lanz = Lanz_abs$Var1, Rivera = Rivera_abs$Var1, Akkouh = Akkouh_abs$Var1,
					FVA_ST = FVA_ST_abs$rxnList, MTA_ST = MTA_ST_abs$del_rxnID_KO,
					
					Lanz = Lanz_norm_t1$Var1, Rivera = Rivera_norm_t1$Var1, Akkouh = Akkouh_norm_t1$Var1,
					FVA_ST = FVA_ST_norm_t1$rxnList, MTA_ST = MTA_ST_norm_t1$del_rxnID_KO,
					
					Lanz = Lanz_norm_t2$Var1, Rivera = Rivera_norm_t2$Var1, Akkouh = Akkouh_norm_t2$Var1,
					FVA_ST = FVA_ST_norm_t2$rxnList, MTA_ST = MTA_ST_norm_t2$del_rxnID_KO)

		y = list(All_Primary_abs$Rxn, All_Ctrl_abs$Rxn, All_HT_abs$Rxn, All_ST_abs$Rxn,
					All_Primary_norm_t1$Rxn, All_Ctrl_norm_t1$Rxn, All_HT_norm_t1$Rxn, All_ST_norm_t1$Rxn,
					All_Primary_norm_t2$Rxn, All_Ctrl_norm_t2$Rxn, All_HT_norm_t2$Rxn, All_ST_norm_t2$Rxn)

		yInt = Reduce(intersect,y); background = length(yInt);
		
#~ 		M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
#~ 		M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
#~ 		rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
#~ 		pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", #display_numbers = TRUE, number_format = "%.1e",
#~ 		main = 'All - ST vs Li+', #color=colorRampPalette(c("green", "yellow", "orangered"))(50))
#~ 		color = c("green", "orange"), breaks = c(0, 0.05, 1), legend = F)
		
#~ 		M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;

###################################
# BD significance across methods
###################################

	# comment this section to plot "# of rxns per subsystem" 
	
		FVA_ST_abs = FVA_ST_abs %>% count(subSystem)
		MTA_ST_abs = MTA_ST_abs %>% count(subSystem_CS)

		FVA_ST_norm_t1 = FVA_ST_norm_t1 %>% count(subSystem)
		MTA_ST_norm_t1 = MTA_ST_norm_t1 %>% count(subSystem_CS)
		
		FVA_ST_norm_t2 = FVA_ST_norm_t2 %>% count(subSystem)
		MTA_ST_norm_t2 = MTA_ST_norm_t2 %>% count(subSystem_CS)

	# model_abs
	
		mat = lst(FVA_ST_abs$subSystem,  MTA_ST_abs$subSystem_CS) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_abs.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_abs.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_ST", "MTA_ST")

		mat = subset(mat, TRUE, c("subSystem", "FVA_ST", "MTA_ST"))

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		pheatmap(mat, cluster_rows=T, cluster_cols=F, display_numbers = round(mat,2), main = 'Abs', 
					color = c("white", "lightblue", "orange"), breaks = c(0, 0.99, 50, max(mat)), legend = F)

	# model_norm_t1
	
		mat = lst(FVA_ST_norm_t1$subSystem,  MTA_ST_norm_t1$subSystem_CS) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_norm_t1.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_norm_t1.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_ST", "MTA_ST")

		mat = subset(mat, TRUE, c("subSystem", "FVA_ST", "MTA_ST"))

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		pheatmap(mat, cluster_rows=T, cluster_cols=F, display_numbers = round(mat,2), main = 'Norm_t1', 
					color = c("white", "lightblue", "orange"), breaks = c(0, 0.99, 50, max(mat)), legend = F)

	# model_norm_t2
	
		mat = lst(FVA_ST_norm_t2$subSystem,  MTA_ST_norm_t2$subSystem_CS) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_norm_t2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_norm_t2.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_ST", "MTA_ST")

		mat = subset(mat, TRUE, c("subSystem", "FVA_ST", "MTA_ST"))

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		pheatmap(mat, cluster_rows=T, cluster_cols=F, display_numbers = round(mat,2), main = 'Norm_t2', 
					color = c("white", "lightblue", "orange"), breaks = c(0, 0.99, 50, max(mat)), legend = F)

	# model_all [i.e., ST]

		mat = lst(FVA_ST_abs$subSystem,  MTA_ST_abs$subSystem_CS,
						FVA_ST_norm_t1$subSystem,  MTA_ST_norm_t1$subSystem_CS,
						FVA_ST_norm_t2$subSystem,  MTA_ST_norm_t2$subSystem_CS) %>% 

		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_st.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_st.csv", header = T, sep = "\t")
		
		colnames(mat) = c("subSystem", "FVA_ST_abs", "FVA_ST_norm_t1", "FVA_ST_norm_t2", 
									"MTA_ST_abs", "MTA_ST_norm_t1", "MTA_ST_norm_t2")		
		
		mat = subset(mat, TRUE, c("subSystem", "FVA_ST_abs", "MTA_ST_abs", 
													"FVA_ST_norm_t1", "MTA_ST_norm_t1", 
													"FVA_ST_norm_t2", "MTA_ST_norm_t2"))
		
		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		pheatmap(mat, cluster_rows=T, cluster_cols=F, main = 'ST', fontsize_row = 8, fontsize_col = 10,
					color = c("white", "lightblue", "orange"), breaks = c(0, 0.99, 50, max(mat)), legend = F)

		mat = cbind(mat, Total = rowSums(mat!= 0))
		mat = mat[order(mat$Total),]
		dotchart(mat$Total, labels = row.names(mat), cex = 0.7, bg = "blue", 
			xlab = "Number of disrupted modules", main = "ST")

		# filter subSystems disrupted by 2 or more methods
		keep = rownames(mat)[rowSums(mat)>2];
		mat_keep = mat[(row.names(mat) %in% keep),]
		write.table(mat_keep, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_st_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)	

###########################
# subSystem - Module matrix
###########################

	# tbl_filt
		ST_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/2.Koskuvi/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_st_filt.csv", header = T, sep = "\t")

				# st	
					df = data.frame("subSystem" = ST_tbl$X, 
								"disruptedBy" = c(ST_tbl$Total))
								
					st_lst = rep(df$subSystem, df$disruptedBy)

		mat = lst(st_lst) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)

		mat = t(mat)
		
		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_mat_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_mat_filt.csv", header = T, sep = "\t")	
		
		mat <- mat[, c("name", "st_lst")]
		
		colnames(mat) = c("subSystem", "ST")
		
		mm <- melt(mat, id="subSystem")

		ggplot(mm, aes(x=reorder(subSystem, -value), y=value, fill=variable)) + 
			geom_bar(stat="identity", color="black", width = 0.7) + theme_classic() +
			scale_fill_manual(values=c("#339cff")) + 
			coord_flip() +
			scale_y_continuous(breaks=c(0:1:10)) + theme(aspect.ratio=1) +
			xlab("") + ylab("# of Disrupted Modules") + labs(fill = "Phenotype") +
			theme(axis.text.y=element_text(size=rel(1.1)))	

###############################################
# Backtracking - subSystems to Modules
###############################################
		
	# st	
	
		DF = ST_tbl
		
		DF$DisruptedModule <- simplify2array(apply(DF[2:7], 1, function(x) paste(names(DF[2:7])[x != 0], collapse = " ")))

		DF$Phenotype <- rep(c("ST"),times=nrow(DF))
		
		DF_ST = DF
		
		ST <- DF_ST %>% relocate(Total, DisruptedModule, .before = FVA_ST_abs)
		
		ST <- ST %>% relocate(Phenotype, .before = X)

		setnames(ST, old = c('X','Total', 'DisruptedModule'), new = c('subSystem','No.of.disrupted.modules', 'disrupted.modules'))

###################
# Load rxn data
###################

	# model_abs
		# results_koskuvi
		FVA_ST_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_significant_abs/ST_rxns_fdr_abs.csv", header = T, sep = "\t")
		MTA_ST_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_prctile_top_significant_abs/ST_rxns_fdr.csv", header = T, sep = "\t")
		
	# model_norm_t1
		# results_koskuvi
		FVA_ST_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_significant_norm_t1/ST_rxns_fdr_norm_t1.csv", header = T, sep = "\t")
		MTA_ST_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_prctile_top_significant_norm_t1/ST_rxns_fdr.csv", header = T, sep = "\t")

	# model_norm_t2
		# results_koskuvi
		FVA_ST_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/FSr_ST/PlotResults/st_tbl_significant_norm_t2/ST_rxns_fdr_norm_t2.csv", header = T, sep = "\t")
		MTA_ST_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_prctile_top_significant_norm_t2/ST_rxns_fdr.csv", header = T, sep = "\t")

###############################################
# Backtracking - Modules to Rxns
###############################################

		ST[,c(1:4)]

		x = strsplit(ST$disrupted.modules, "\\s+")
		
	# Nucleotide interconversion (ST_S1)
	
		DF_d1 = get(x[[1]][1])
		DF_d2 = get(x[[1]][2])
		DF_d3 = get(x[[1]][3])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem_CS == toString(ST$subSystem[[1]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CS == toString(ST$subSystem[[1]]))
		DF_d3 = DF_d3 %>% dplyr::filter(subSystem_CS == toString(ST$subSystem[[1]]))
		DF_d1$disrupted.module <- rep(c(x[[1]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[1]][2]),times=nrow(DF_d2))
		DF_d3$disrupted.module <- rep(c(x[[1]][3]),times=nrow(DF_d3))

		DF_d1 <- subset(DF_d1, TRUE, c("subSystem_CS", "disrupted.module", "del_rxnID_KO"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CS", "disrupted.module", "del_rxnID_KO"))
		DF_d3 <- subset(DF_d3, TRUE, c("subSystem_CS", "disrupted.module", "del_rxnID_KO"))
		colnames(DF_d1) = c("subSystem", "disrupted.module", "rxnList")
		colnames(DF_d2) = c("subSystem", "disrupted.module", "rxnList")
		colnames(DF_d3) = c("subSystem", "disrupted.module", "rxnList")
		ST_S1 = RBIND(list(DF_d1, DF_d2, DF_d3))

		write.table(ST_S1, "PlotResults/plotDataOverlaysHyper_Tbl_Final/ST_S1.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

	# Fatty acid oxidation (ST_S2)
	
		DF_d1 = get(x[[2]][1])
		DF_d2 = get(x[[2]][2])
		DF_d3 = get(x[[2]][3])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem == toString(ST$subSystem[[2]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem == toString(ST$subSystem[[2]]))
		DF_d3 = DF_d3 %>% dplyr::filter(subSystem_CS == toString(ST$subSystem[[2]]))
		DF_d1$disrupted.module <- rep(c(x[[2]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[2]][2]),times=nrow(DF_d2))
		DF_d3$disrupted.module <- rep(c(x[[2]][3]),times=nrow(DF_d3))

		DF_d1 <- subset(DF_d1, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d3 <- subset(DF_d3, TRUE, c("subSystem_CS", "disrupted.module", "del_rxnID_KO"))
		colnames(DF_d3) = c("subSystem", "disrupted.module", "rxnList")
		ST_S2 = RBIND(list(DF_d1, DF_d2, DF_d3))
		
		write.table(ST_S2, "PlotResults/plotDataOverlaysHyper_Tbl_Final/ST_S2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)		
		
###############################################
# Backtracking - Minerva Visualization
###############################################

######################
# Final-subSystem-wise (n=2)
######################

			# duplicate ST files
			DF_1 = ST_S1
			DF_2 = ST_S2

			# add 'r_' as prefix to Rxns
			DF_1$rxnList = paste0('r_', DF_1$rxnList)
			DF_2$rxnList = paste0('r_', DF_2$rxnList)
		
			# rm duplicates
			DF_1 <- DF_1[!duplicated(DF_1[,c('rxnList')]),]
			DF_2 <- DF_2[!duplicated(DF_2[,c('rxnList')]),]	
			
			# add 'lineWidth'
			DF_1$lineWidth <- rep(c("3"),times=nrow(DF_1))
			DF_2$lineWidth <- rep(c("3"),times=nrow(DF_2))			

			# add 'color'
			DF_1$color <- rep(c("#0000FF"),times=nrow(DF_1))	#blue
			DF_2$color <- rep(c("#0000FF"),times=nrow(DF_2))	#blue

			# slice data
			DF_1 <- subset(DF_1, TRUE, c("rxnList", "lineWidth", "color"))
			DF_2 <- subset(DF_2, TRUE, c("rxnList", "lineWidth", "color"))

			# rename colnames
			colnames(DF_1) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_2) = c("reactionIdentifier", "lineWidth", "color")

			# write table
			write.table(DF_1, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_1_ST_S1_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			write.table(DF_2, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_2_ST_S2_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			
######################
# Final-phenotype-wise (n=1)
######################		
		
			# Merge DFs for Minerva
			DF1 = RBIND(list(ST_S1, ST_S2))			
				
			# add 'r_' as prefix to Rxns
			DF1$rxnList = paste0('r_', DF1$rxnList)
		
			# rm duplicates
			DF1 <- DF1[!duplicated(DF1[,c('rxnList')]),]

			# add 'lineWidth'
			DF1$lineWidth <- rep(c("3"),times=nrow(DF1))

			# add 'color'
			DF1$color <- rep(c("#0000FF"),times=nrow(DF1))	#blue

			# slice data
			DF1 <- subset(DF1, TRUE, c("rxnList", "lineWidth", "color"))

			# rename colnames
			colnames(DF1) = c("reactionIdentifier", "lineWidth", "color")

			dim(DF1)

			# write table
			write.table(DF1, "PlotResults/ST_Rxns.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)







