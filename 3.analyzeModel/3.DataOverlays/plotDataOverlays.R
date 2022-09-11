#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table) 

setwd("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays")

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

###################
# Load rxn data
###################

	# model_abs
		# background
		All_Primary_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_Primary_TP_abs.csv", header = T, sep = "\t")
		All_Ctrl_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_Ctrl_TP_abs.csv", header = T, sep = "\t")
		All_BD_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_TP_abs.csv", header = T, sep = "\t")
		All_BD_R_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_R_TP_abs.csv", header = T, sep = "\t")
		All_BD_NR_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_NR_TP_abs.csv", header = T, sep = "\t")
		# results_vadodaria
		FVA_BD_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_abs/bd_lumped_rxns_fdr_abs.csv", header = T, sep = "\t")
		FVA_BD_R_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_abs/bd_r_rxns_fdr_abs.csv", header = T, sep = "\t")
		FVA_BD_NR_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_abs/bd_nr_rxns_fdr_abs.csv", header = T, sep = "\t")
		MTA_BD_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_abs/bd_lumped_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_R_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_abs/bd_r_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_NR_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_abs/bd_nr_rxns_fdr.csv", header = T, sep = "\t")
		# results_lithium
		Lanz_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/AnalysisFiles/Lanz_vs_Ctrl_abs.csv", header = T, sep = "\t")
		Rivera_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/AnalysisFiles/Rivera_vs_Ctrl_abs.csv", header = T, sep = "\t")
		
	# model_norm_t1
		# background
		All_Primary_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_Primary_TP_norm_t1.csv", header = T, sep = "\t")
		All_Ctrl_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_Ctrl_TP_norm_t1.csv", header = T, sep = "\t")
		All_BD_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_TP_norm_t1.csv", header = T, sep = "\t")
		All_BD_R_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_R_TP_norm_t1.csv", header = T, sep = "\t")
		All_BD_NR_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_NR_TP_norm_t1.csv", header = T, sep = "\t")
		# results_vadodaria
		FVA_BD_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t1/bd_lumped_rxns_fdr_norm_t1.csv", header = T, sep = "\t")
		FVA_BD_R_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t1/bd_r_rxns_fdr_norm_t1.csv", header = T, sep = "\t")
		FVA_BD_NR_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t1/bd_nr_rxns_fdr_norm_t1.csv", header = T, sep = "\t")
		MTA_BD_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t1/bd_lumped_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_R_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t1/bd_r_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_NR_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t1/bd_nr_rxns_fdr.csv", header = T, sep = "\t")
		# results_lithium
		Lanz_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/AnalysisFiles/Lanz_vs_Ctrl_norm_t1.csv", header = T, sep = "\t")
		Rivera_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/AnalysisFiles/Rivera_vs_Ctrl_norm_t1.csv", header = T, sep = "\t")

	# model_norm_t2
		# background
		All_Primary_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_Primary_TP_norm_t2.csv", header = T, sep = "\t")
		All_Ctrl_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_Ctrl_TP_norm_t2.csv", header = T, sep = "\t")
		All_BD_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_TP_norm_t2.csv", header = T, sep = "\t")
		All_BD_R_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_R_TP_norm_t2.csv", header = T, sep = "\t")
		All_BD_NR_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_NR_TP_norm_t2.csv", header = T, sep = "\t")
		# results_vadodaria
		FVA_BD_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t2/bd_lumped_rxns_fdr_norm_t2.csv", header = T, sep = "\t")
		FVA_BD_R_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t2/bd_r_rxns_fdr_norm_t2.csv", header = T, sep = "\t")
		FVA_BD_NR_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t2/bd_nr_rxns_fdr_norm_t2.csv", header = T, sep = "\t")
		MTA_BD_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t2/bd_lumped_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_R_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t2/bd_r_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_NR_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t2/bd_nr_rxns_fdr.csv", header = T, sep = "\t")
		# results_lithium
		Lanz_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/AnalysisFiles/Lanz_vs_Ctrl_norm_t2.csv", header = T, sep = "\t")
		Rivera_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/AnalysisFiles/Rivera_vs_Ctrl_norm_t2.csv", header = T, sep = "\t")

############################
# BD vs Lithium
############################
	
	# model_abs

		x = list(Lanz = Lanz_abs$Var1, Rivera = Rivera_abs$Var1,
					FVA_BD = FVA_BD_abs$rxnList, FVA_BD_R = FVA_BD_R_abs$rxnList, FVA_BD_NR = FVA_BD_NR_abs$rxnList,
						MTA_BD = MTA_BD_abs$del_rxnID_KO, MTA_BD_R = MTA_BD_R_abs$del_rxnID_KO, MTA_BD_NR = MTA_BD_NR_abs$del_rxnID_KO)

		y = list(All_Primary_abs$Rxn, All_Ctrl_abs$Rxn, All_BD_abs$Rxn, All_BD_R_abs$Rxn, All_BD_NR_abs$Rxn)
		yInt = Reduce(intersect,y); background = length(yInt);
		
		M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
		M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
		rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
		pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
		main = 'Abs - BD vs Li+', #color=colorRampPalette(c("green", "yellow", "orangered"))(50))
		color = c("green", "orange"), breaks = c(0, 0.05, 1), legend = F)

		M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
#~ 		M.int
	
	# model_norm_t1

		x = list(Lanz = Lanz_norm_t1$Var1, Rivera = Rivera_norm_t1$Var1,
					FVA_BD = FVA_BD_norm_t1$rxnList, FVA_BD_R = FVA_BD_R_norm_t1$rxnList, FVA_BD_NR = FVA_BD_NR_norm_t1$rxnList,
						MTA_BD = MTA_BD_norm_t1$del_rxnID_KO, MTA_BD_R = MTA_BD_R_norm_t1$del_rxnID_KO, MTA_BD_NR = MTA_BD_NR_norm_t1$del_rxnID_KO)

		y = list(All_Primary_norm_t1$Rxn, All_Ctrl_norm_t1$Rxn, All_BD_norm_t1$Rxn, All_BD_R_norm_t1$Rxn, All_BD_NR_norm_t1$Rxn)
		yInt = Reduce(intersect,y); background = length(yInt);
		
		M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
		M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
		rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
		pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
		main = 'Norm_T1 - BD vs Li+', #color=colorRampPalette(c("green", "yellow", "orangered"))(50))
		color = c("green", "orange"), breaks = c(0, 0.05, 1), legend = F)
		
		M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
#~ 		M.int
	
	# model_norm_t2

		x = list(Lanz = Lanz_norm_t2$Var1, Rivera = Rivera_norm_t2$Var1,
					FVA_BD = FVA_BD_norm_t2$rxnList, FVA_BD_R = FVA_BD_R_norm_t2$rxnList, FVA_BD_NR = FVA_BD_NR_norm_t2$rxnList,
						MTA_BD = MTA_BD_norm_t2$del_rxnID_KO, MTA_BD_R = MTA_BD_R_norm_t2$del_rxnID_KO, MTA_BD_NR = MTA_BD_NR_norm_t2$del_rxnID_KO)

		y = list(All_Primary_norm_t2$Rxn, All_Ctrl_norm_t2$Rxn, All_BD_norm_t2$Rxn, All_BD_R_norm_t2$Rxn, All_BD_NR_norm_t2$Rxn)
		yInt = Reduce(intersect,y); background = length(yInt);
		
		M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
		M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
		rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
		pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
		main = 'Norm_T2 - BD vs Li+', #color=colorRampPalette(c("green", "yellow", "orangered"))(50))
		color = c("green", "orange"), breaks = c(0, 0.05, 1), legend = F)
		
		M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
#~ 		M.int

	# model_all

		x = list(Lanz_abs = Lanz_abs$Var1, Rivera_abs = Rivera_abs$Var1,
					FVA_BD_abs = FVA_BD_abs$rxnList, FVA_BD_R_abs = FVA_BD_R_abs$rxnList, FVA_BD_NR_abs = FVA_BD_NR_abs$rxnList,
						MTA_BD_abs = MTA_BD_abs$del_rxnID_KO, MTA_BD_R_abs = MTA_BD_R_abs$del_rxnID_KO, MTA_BD_NR_abs = MTA_BD_NR_abs$del_rxnID_KO,
				
				Lanz_norm_t1 = Lanz_norm_t1$Var1, Rivera_norm_t1 = Rivera_norm_t1$Var1,
					FVA_BD_norm_t1 = FVA_BD_norm_t1$rxnList, FVA_BD_R_norm_t1 = FVA_BD_R_norm_t1$rxnList, FVA_BD_NR_norm_t1 = FVA_BD_NR_norm_t1$rxnList,
						MTA_BD_norm_t1 = MTA_BD_norm_t1$del_rxnID_KO, MTA_BD_R_norm_t1 = MTA_BD_R_norm_t1$del_rxnID_KO, MTA_BD_NR_norm_t1 = MTA_BD_NR_norm_t1$del_rxnID_KO,
						
				Lanz_norm_t2 = Lanz_norm_t2$Var1, Rivera_norm_t2 = Rivera_norm_t2$Var1,
					FVA_BD_norm_t2 = FVA_BD_norm_t2$rxnList, FVA_BD_R_norm_t2 = FVA_BD_R_norm_t2$rxnList, FVA_BD_NR_norm_t2 = FVA_BD_NR_norm_t2$rxnList,
						MTA_BD_norm_t2 = MTA_BD_norm_t2$del_rxnID_KO, MTA_BD_R_norm_t2 = MTA_BD_R_norm_t2$del_rxnID_KO, MTA_BD_NR_norm_t2 = MTA_BD_NR_norm_t2$del_rxnID_KO)

		y = list(All_Primary_abs$Rxn, All_Ctrl_abs$Rxn, All_BD_abs$Rxn, All_BD_R_abs$Rxn, All_BD_NR_abs$Rxn,
					All_Primary_norm_t1$Rxn, All_Ctrl_norm_t1$Rxn, All_BD_norm_t1$Rxn, All_BD_R_norm_t1$Rxn, All_BD_NR_norm_t1$Rxn,
						All_Primary_norm_t2$Rxn, All_Ctrl_norm_t2$Rxn, All_BD_norm_t2$Rxn, All_BD_R_norm_t2$Rxn, All_BD_NR_norm_t2$Rxn)
						
		yInt = Reduce(intersect,y); background = length(yInt);
		
		M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
		M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
		rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
		pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", #display_numbers = TRUE, number_format = "%.1e",
		main = 'All - BD vs Li+', #color=colorRampPalette(c("green", "yellow", "orangered"))(50))
		color = c("green", "orange"), breaks = c(0, 0.05, 1), legend = F)
		
		M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
#~ 		M.int

###################################
# BD significance across methods
###################################

	# comment this section to plot "# of rxns per subsystem" 
	
		FVA_BD_abs = FVA_BD_abs %>% count(subSystem)
		FVA_BD_R_abs = FVA_BD_R_abs %>% count(subSystem)
		FVA_BD_NR_abs = FVA_BD_NR_abs %>% count(subSystem)
		MTA_BD_abs = MTA_BD_abs %>% count(subSystem_CB)
		MTA_BD_R_abs = MTA_BD_R_abs %>% count(subSystem_CBR)
		MTA_BD_NR_abs = MTA_BD_NR_abs %>% count(subSystem_CBNR)
		
		FVA_BD_norm_t1 = FVA_BD_norm_t1 %>% count(subSystem)
		FVA_BD_R_norm_t1 = FVA_BD_R_norm_t1 %>% count(subSystem)
		FVA_BD_NR_norm_t1 = FVA_BD_NR_norm_t1 %>% count(subSystem)
		MTA_BD_norm_t1 = MTA_BD_norm_t1 %>% count(subSystem_CB)
		MTA_BD_R_norm_t1 = MTA_BD_R_norm_t1 %>% count(subSystem_CBR)
		MTA_BD_NR_norm_t1 = MTA_BD_NR_norm_t1 %>% count(subSystem_CBNR)
		
		FVA_BD_norm_t2 = FVA_BD_norm_t2 %>% count(subSystem)
		FVA_BD_R_norm_t2 = FVA_BD_R_norm_t2 %>% count(subSystem)
		FVA_BD_NR_norm_t2 = FVA_BD_NR_norm_t2 %>% count(subSystem)
		MTA_BD_norm_t2 = MTA_BD_norm_t2 %>% count(subSystem_CB)
		MTA_BD_R_norm_t2 = MTA_BD_R_norm_t2 %>% count(subSystem_CBR)
		MTA_BD_NR_norm_t2 = MTA_BD_NR_norm_t2 %>% count(subSystem_CBNR)

	# model_abs
	
		mat = lst(FVA_BD_abs$subSystem, FVA_BD_R_abs$subSystem, FVA_BD_NR_abs$subSystem, 
					MTA_BD_abs$subSystem_CB, MTA_BD_R_abs$subSystem_CBR, MTA_BD_NR_abs$subSystem_CBNR) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_abs.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_abs.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_BD", "FVA_BD_R", "FVA_BD_NR", 
							"MTA_BD", "MTA_BD_R", "MTA_BD_NR")

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		pheatmap(mat, cluster_rows=T, cluster_cols=F, display_numbers = round(mat,2), main = 'Abs', 
					color = c("white", "lightblue", "orange"), breaks = c(0, 0.99, 50, max(mat)), legend = F)

	# model_norm_t1
	
		mat = lst(FVA_BD_norm_t1$subSystem, FVA_BD_R_norm_t1$subSystem, FVA_BD_NR_norm_t1$subSystem, 
					MTA_BD_norm_t1$subSystem_CB, MTA_BD_R_norm_t1$subSystem_CBR, MTA_BD_NR_norm_t1$subSystem_CBNR) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_norm_t1.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_norm_t1.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_BD", "FVA_BD_R", "FVA_BD_NR", 
							"MTA_BD", "MTA_BD_R", "MTA_BD_NR")

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		pheatmap(mat, cluster_rows=T, cluster_cols=F, display_numbers = round(mat,2), main = 'Norm_T1', 
					color = c("white", "lightblue", "orange"), breaks = c(0, 0.99, 50, max(mat)), legend = F)

	# model_norm_t2
	
		mat = lst(FVA_BD_norm_t2$subSystem, FVA_BD_R_norm_t2$subSystem, FVA_BD_NR_norm_t2$subSystem, 
					MTA_BD_norm_t2$subSystem_CB, MTA_BD_R_norm_t2$subSystem_CBR, MTA_BD_NR_norm_t2$subSystem_CBNR) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_norm_t2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_norm_t2.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_BD", "FVA_BD_R", "FVA_BD_NR", 
							"MTA_BD", "MTA_BD_R", "MTA_BD_NR")

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		pheatmap(mat, cluster_rows=T, cluster_cols=F, display_numbers = round(mat,2), main = 'Norm_T2', 
					color = c("white", "lightblue", "orange"), breaks = c(0, 0.99, 50, max(mat)), legend = F)
					

	# model_all

		mat = lst(FVA_BD_abs$subSystem, FVA_BD_R_abs$subSystem, FVA_BD_NR_abs$subSystem, 
					MTA_BD_abs$subSystem_CB, MTA_BD_R_abs$subSystem_CBR, MTA_BD_NR_abs$subSystem_CBNR,
					FVA_BD_norm_t1$subSystem, FVA_BD_R_norm_t1$subSystem, FVA_BD_NR_norm_t1$subSystem, 
					MTA_BD_norm_t1$subSystem_CB, MTA_BD_R_norm_t1$subSystem_CBR, MTA_BD_NR_norm_t1$subSystem_CBNR,
					FVA_BD_norm_t2$subSystem, FVA_BD_R_norm_t2$subSystem, FVA_BD_NR_norm_t2$subSystem, 
					MTA_BD_norm_t2$subSystem_CB, MTA_BD_R_norm_t2$subSystem_CBR, MTA_BD_NR_norm_t2$subSystem_CBNR) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_all.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_all.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_BD_abs", "FVA_BD_R_abs", "FVA_BD_NR_abs", "MTA_BD_abs", "MTA_BD_R_abs", "MTA_BD_NR_abs",
							"FVA_BD_norm_t1", "FVA_BD_R_norm_t1", "FVA_BD_NR_norm_t1", "MTA_BD_norm_t1", "MTA_BD_R_norm_t1", "MTA_BD_NR_norm_t1",
							"FVA_BD_norm_t2", "FVA_BD_R_norm_t2", "FVA_BD_NR_norm_t2", "MTA_BD_norm_t2", "MTA_BD_R_norm_t2", "MTA_BD_NR_norm_t2")

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		pheatmap(mat, cluster_rows=T, cluster_cols=F, main = 'All', fontsize_row = 8, fontsize_col = 10,
					color = c("white", "lightblue", "orange"), breaks = c(0, 0.99, 50, max(mat)), legend = F)

		mat = cbind(mat, Total = rowSums(mat!= 0))
		mat = mat[order(mat$Total),]
		dotchart(mat$Total, labels = row.names(mat), cex = 0.7, bg = "blue", 
			xlab = "Number of disrupted modules", main = "All")
			

	# model_bd

		mat = lst(FVA_BD_abs$subSystem,
					MTA_BD_abs$subSystem_CB,
					FVA_BD_norm_t1$subSystem,
					MTA_BD_norm_t1$subSystem_CB,
					FVA_BD_norm_t2$subSystem,
					MTA_BD_norm_t2$subSystem_CB) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_BD_abs", "FVA_BD_norm_t1", 
							"FVA_BD_norm_t2", "MTA_BD_abs", 
							"MTA_BD_norm_t1", "MTA_BD_norm_t2")

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		mat = cbind(mat, Total = rowSums(mat!= 0))
		mat = mat[order(mat$Total),]
		dotchart(mat$Total, labels = row.names(mat), cex = 0.7, bg = "blue", 
			xlab = "Number of disrupted modules", main = "BD_Lumped")		
			
		# filter subSystems disrupted by 2 or more methods
		keep = rownames(mat)[rowSums(mat)>2];
		mat_keep = mat[(row.names(mat) %in% keep),]
		write.table(mat_keep, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)	
		
	# model_bd_r

		mat = lst(FVA_BD_R_abs$subSystem,
					MTA_BD_R_abs$subSystem_CBR,
					FVA_BD_R_norm_t1$subSystem,
					MTA_BD_R_norm_t1$subSystem_CBR,
					FVA_BD_R_norm_t2$subSystem,
					MTA_BD_R_norm_t2$subSystem_CBR) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_r.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_r.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_BD_R_abs", "FVA_BD_R_norm_t1", 
							"FVA_BD_R_norm_t2", "MTA_BD_R_abs", 
							"MTA_BD_R_norm_t1", "MTA_BD_R_norm_t2")

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		mat = cbind(mat, Total = rowSums(mat!= 0))
		mat = mat[order(mat$Total),]
		dotchart(mat$Total, labels = row.names(mat), cex = 0.7, bg = "blue", 
			xlab = "Number of disrupted modules", main = "BD_R")
			
		# filter subSystems disrupted by 2 or more methods
		keep = rownames(mat)[rowSums(mat)>2];
		mat_keep = mat[(row.names(mat) %in% keep),]
		write.table(mat_keep, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_r_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)	

	# model_bd_nr

		mat = lst(FVA_BD_NR_abs$subSystem,
					MTA_BD_NR_abs$subSystem_CBNR,
					FVA_BD_NR_norm_t1$subSystem,
					MTA_BD_NR_norm_t1$subSystem_CBNR,
					FVA_BD_NR_norm_t2$subSystem,
					MTA_BD_NR_norm_t2$subSystem_CBNR) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)

		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_nr.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_nr.csv", header = T, sep = "\t")

		colnames(mat) = c("subSystem", "FVA_BD_NR_abs", "FVA_BD_NR_norm_t1", 
							"FVA_BD_NR_norm_t2", "MTA_BD_NR_abs", 
							"MTA_BD_NR_norm_t1", "MTA_BD_NR_norm_t2")

		row.names(mat) <- mat$subSystem
		
		mat = mat[,colnames(mat)!="subSystem"]

		mat = cbind(mat, Total = rowSums(mat!= 0))
		mat = mat[order(mat$Total),]
		dotchart(mat$Total, labels = row.names(mat), cex = 0.7, bg = "blue",
			xlab = "Number of disrupted modules", main = "BD_NR")
			
		# filter subSystems disrupted by 2 or more methods
		keep = rownames(mat)[rowSums(mat)>2];
		mat_keep = mat[(row.names(mat) %in% keep),]
		write.table(mat_keep, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_nr_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
		
		
###########################
# subSystem - Module matrix
###########################

	# tbl_filt
		BD_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_filt.csv", header = T, sep = "\t")
		BD_R_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_r_filt.csv", header = T, sep = "\t")
		BD_NR_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/plotDataOverlaysHyper_Tbl/Tbl_bd_nr_filt.csv", header = T, sep = "\t")
		
				# bd	
					df = data.frame("subSystem" = BD_tbl$X, 
								"disruptedBy" = c(BD_tbl$Total))
								
					bd_lst = rep(df$subSystem, df$disruptedBy)


				# bd_r
					df = data.frame("subSystem" = BD_R_tbl$X, 
								"disruptedBy" = c(BD_R_tbl$Total))
								
					bd_r_lst = rep(df$subSystem, df$disruptedBy)


				# bd_nr
					df = data.frame("subSystem" = BD_NR_tbl$X, 
								"disruptedBy" = c(BD_NR_tbl$Total))
								
					bd_nr_lst = rep(df$subSystem, df$disruptedBy)
		

		mat = lst(bd_lst,
					bd_r_lst,
					bd_nr_lst) %>% 
					
		  enframe %>% 
		  unnest %>% 
		  count(name, value) %>% 
		  spread(value, n, fill = 0)
		  
		mat = t(mat)
		
		write.table(mat, "PlotResults/plotDataOverlaysHyper_Tbl/Tbl_mat_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/plotDataOverlaysHyper_Tbl/Tbl_mat_filt.csv", header = T, sep = "\t")	
		
		mat <- mat[, c("name", "bd_lst", "bd_r_lst", "bd_nr_lst")]
		
		colnames(mat) = c("subSystem", "BD", "BD_R", "BD_NR")
		
		mm <- melt(mat, id="subSystem")
			
		ggplot(mm, aes(x=reorder(subSystem, -value), y=value, fill=variable)) + 
			geom_bar(stat="identity", color="black", width = 0.7) + theme_classic() +
			scale_fill_manual(values=c("#339cff", "#33ff58", "#ff4233")) + 
			coord_flip() +
			scale_y_continuous(breaks=c(0:1:10)) + theme(aspect.ratio=1) +
			xlab("") + ylab("# of Disrupted Modules") + labs(fill = "Phenotype") +
			theme(axis.text.y=element_text(size=rel(1.1)))	

###############################################
# Backtracking - subSystems to Modules
###############################################
		
	# bd	
	
		DF = BD_tbl
		
		DF$DisruptedModule <- simplify2array(apply(DF[2:7], 1, function(x) paste(names(DF[2:7])[x != 0], collapse = " ")))
		
		DF$Phenotype <- rep(c("BD"),times=nrow(DF))
		
		DF_BD = DF
		
		DF_BD
					
	# bd_r
	
		DF = BD_R_tbl
		
		DF$DisruptedModule <- simplify2array(apply(DF[2:7], 1, function(x) paste(names(DF[2:7])[x != 0], collapse = " ")))
		
		DF$Phenotype <- rep(c("BD_R"),times=nrow(DF))
		
		DF_BD_R = DF

	# bd_nr
	
		DF = BD_NR_tbl
		
		DF$DisruptedModule <- simplify2array(apply(DF[2:7], 1, function(x) paste(names(DF[2:7])[x != 0], collapse = " ")))
		
		DF$Phenotype <- rep(c("BD_NR"),times=nrow(DF))
		
		DF_BD_NR = DF 


	# merge bd, bd_r, bd_nr

		df_list <- list(DF_BD, DF_BD_R, DF_BD_NR)
		
		RBIND <- function(datalist) {
					  require(plyr)
					  temp <- rbind.fill(datalist)
					  temp
					}
					
		BD = RBIND(df_list)
		
		BD <- BD %>% relocate(Total, DisruptedModule, .before = FVA_BD_abs)
		
		BD <- BD %>% relocate(Phenotype, .before = X)

		setnames(BD, old = c('X','Total', 'DisruptedModule'), new = c('subSystem','No.of.disrupted.modules', 'disrupted.modules'))

#~ 		BD

###################
# Load rxn data
###################

	# model_abs
		# results_vadodaria
		FVA_BD_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_abs/bd_lumped_rxns_fdr_abs.csv", header = T, sep = "\t")
		FVA_BD_R_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_abs/bd_r_rxns_fdr_abs.csv", header = T, sep = "\t")
		FVA_BD_NR_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_abs/bd_nr_rxns_fdr_abs.csv", header = T, sep = "\t")
		MTA_BD_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_abs/bd_lumped_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_R_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_abs/bd_r_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_NR_abs <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_abs/bd_nr_rxns_fdr.csv", header = T, sep = "\t")
		
	# model_norm_t1
		# results_vadodaria
		FVA_BD_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t1/bd_lumped_rxns_fdr_norm_t1.csv", header = T, sep = "\t")
		FVA_BD_R_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t1/bd_r_rxns_fdr_norm_t1.csv", header = T, sep = "\t")
		FVA_BD_NR_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t1/bd_nr_rxns_fdr_norm_t1.csv", header = T, sep = "\t")
		MTA_BD_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t1/bd_lumped_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_R_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t1/bd_r_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_NR_norm_t1 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t1/bd_nr_rxns_fdr.csv", header = T, sep = "\t")

	# model_norm_t2
		# results_vadodaria
		FVA_BD_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t2/bd_lumped_rxns_fdr_norm_t2.csv", header = T, sep = "\t")
		FVA_BD_R_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t2/bd_r_rxns_fdr_norm_t2.csv", header = T, sep = "\t")
		FVA_BD_NR_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_significant_norm_t2/bd_nr_rxns_fdr_norm_t2.csv", header = T, sep = "\t")
		MTA_BD_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t2/bd_lumped_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_R_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t2/bd_r_rxns_fdr.csv", header = T, sep = "\t")
		MTA_BD_NR_norm_t2 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_significant_norm_t2/bd_nr_rxns_fdr.csv", header = T, sep = "\t")

###############################################
# Backtracking - Modules to Rxns
###############################################
		
		BD[,c(1:4)]
		
		x = strsplit(BD$disrupted.modules, "\\s+")

	# Methionine and cysteine metabolism (BD_S1)
	
		DF_d1 = get(x[[1]][1])
		DF_d2 = get(x[[1]][2])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem == toString(BD$subSystem[[1]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem == toString(BD$subSystem[[1]]))
		DF_d1$disrupted.module <- rep(c(x[[1]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[1]][2]),times=nrow(DF_d2))
		
		DF = RBIND(list(DF_d1, DF_d2))
		DF <- DF %>% relocate(disrupted.module, .before = rxnList)
		BD_S1 <- DF %>% relocate(subSystem, .before = disrupted.module)
		BD_S1 <- subset(BD_S1, TRUE, c("subSystem", "disrupted.module", "rxnList"))
					
		write.table(BD_S1, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S1.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)


	# Keratan sulfate degradation (BD_S2)
	
		DF_d1 = get(x[[2]][1])
		DF_d2 = get(x[[2]][2])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem == toString(BD$subSystem[[2]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CB == toString(BD$subSystem[[2]]))
		DF_d1$disrupted.module <- rep(c(x[[2]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[2]][2]),times=nrow(DF_d2))
		
		DF_d1 <- subset(DF_d1, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		colnames(DF_d2) = c("subSystem", "disrupted.module", "rxnList")
		BD_S2 = RBIND(list(DF_d1, DF_d2))
				
		write.table(BD_S2, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
		
	# Fatty acid oxidation (BD_S3)
	
		DF_d1 = get(x[[3]][1])
		DF_d2 = get(x[[3]][2])
		DF_d3 = get(x[[3]][3])
		DF_d4 = get(x[[3]][4])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem == toString(BD$subSystem[[3]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem == toString(BD$subSystem[[3]]))
		DF_d3 = DF_d3 %>% dplyr::filter(subSystem_CB == toString(BD$subSystem[[3]]))
		DF_d4 = DF_d4 %>% dplyr::filter(subSystem_CB == toString(BD$subSystem[[3]]))
		DF_d1$disrupted.module <- rep(c(x[[3]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[3]][2]),times=nrow(DF_d2))
		DF_d3$disrupted.module <- rep(c(x[[3]][3]),times=nrow(DF_d3))
		DF_d4$disrupted.module <- rep(c(x[[3]][4]),times=nrow(DF_d4))

		DF_d1 <- subset(DF_d1, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d3 <- subset(DF_d3, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		DF_d4 <- subset(DF_d4, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		colnames(DF_d3) = c("subSystem", "disrupted.module", "rxnList")
		colnames(DF_d4) = c("subSystem", "disrupted.module", "rxnList")
		BD_S3 = RBIND(list(DF_d1, DF_d2, DF_d3, DF_d4))
		
		write.table(BD_S3, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S3.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
		
	# ROS detoxification (BD_R_S1)
	
		DF_d1 = get(x[[4]][1])
		DF_d2 = get(x[[4]][2])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem == toString(BD$subSystem[[4]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CBR == toString(BD$subSystem[[4]]))
		DF_d1$disrupted.module <- rep(c(x[[4]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[4]][2]),times=nrow(DF_d2))
		
		DF_d1 <- subset(DF_d1, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))		
		colnames(DF_d2) = c("subSystem", "disrupted.module", "rxnList")
		BD_R_S1 = RBIND(list(DF_d1, DF_d2))
		
		write.table(BD_R_S1, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S1.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
		
	# Drug metabolism (BD_R_S2)
	
		DF_d1 = get(x[[5]][1])
		DF_d2 = get(x[[5]][2])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem_CBR == toString(BD$subSystem[[5]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CBR == toString(BD$subSystem[[5]]))
		DF_d1$disrupted.module <- rep(c(x[[5]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[5]][2]),times=nrow(DF_d2))
		
		DF_d1 <- subset(DF_d1, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		colnames(DF_d1) = c("subSystem", "disrupted.module", "rxnList")
		colnames(DF_d1) = c("subSystem", "disrupted.module", "rxnList")
		BD_R_S2 = RBIND(list(DF_d1, DF_d2))
		
		write.table(BD_R_S2, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

	# Transport, extracellular (BD_R_S3)
	
		DF_d1 = get(x[[6]][1])
		DF_d2 = get(x[[6]][2])
		DF_d3 = get(x[[6]][3])		
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem_CBR == toString(BD$subSystem[[6]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CBR == toString(BD$subSystem[[6]]))
		DF_d3 = DF_d3 %>% dplyr::filter(subSystem_CBR == toString(BD$subSystem[[6]]))
		DF_d1$disrupted.module <- rep(c(x[[6]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[6]][2]),times=nrow(DF_d2))
		DF_d3$disrupted.module <- rep(c(x[[6]][3]),times=nrow(DF_d3))

		DF_d1 <- subset(DF_d1, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		DF_d3 <- subset(DF_d3, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		colnames(DF_d1) = c("subSystem", "disrupted.module", "rxnList")
		colnames(DF_d2) = c("subSystem", "disrupted.module", "rxnList")
		colnames(DF_d3) = c("subSystem", "disrupted.module", "rxnList")
		BD_R_S3 = RBIND(list(DF_d1, DF_d2, DF_d3))

		write.table(BD_R_S3, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S3.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

	# Cholesterol metabolism (BD_NR_S1)
	
		DF_d1 = get(x[[7]][1])
		DF_d2 = get(x[[7]][2])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem == toString(BD$subSystem[[7]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CBNR == toString(BD$subSystem[[7]]))
		DF_d1$disrupted.module <- rep(c(x[[7]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[7]][2]),times=nrow(DF_d2))	
				
		DF_d1 <- subset(DF_d1, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))	
		colnames(DF_d2) = c("subSystem", "disrupted.module", "rxnList")
		BD_NR_S1 = RBIND(list(DF_d1, DF_d2))

		write.table(BD_NR_S1, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S1.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

	# Inositol phosphate metabolism (BD_NR_S2)
	
		DF_d1 = get(x[[8]][1])
		DF_d2 = get(x[[8]][2])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem == toString(BD$subSystem[[8]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CBNR == toString(BD$subSystem[[8]]))
		DF_d1$disrupted.module <- rep(c(x[[8]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[8]][2]),times=nrow(DF_d2))	
				
		DF_d1 <- subset(DF_d1, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))	
		colnames(DF_d2) = c("subSystem", "disrupted.module", "rxnList")
		BD_NR_S2 = RBIND(list(DF_d1, DF_d2))

		write.table(BD_NR_S2, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

	# Fatty acid oxidation (BD_NR_S3)
	
		DF_d1 = get(x[[9]][1])
		DF_d2 = get(x[[9]][2])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem == toString(BD$subSystem[[9]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CBNR == toString(BD$subSystem[[9]]))
		DF_d1$disrupted.module <- rep(c(x[[9]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[9]][2]),times=nrow(DF_d2))	
				
		DF_d1 <- subset(DF_d1, TRUE, c("subSystem", "disrupted.module", "rxnList"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))	
		colnames(DF_d2) = c("subSystem", "disrupted.module", "rxnList")
		BD_NR_S3 = RBIND(list(DF_d1, DF_d2))

		write.table(BD_NR_S3, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S3.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

	# Heparan sulfate degradation (BD_NR_S4)
	
		DF_d1 = get(x[[10]][1])
		DF_d2 = get(x[[10]][2])
		DF_d1 = DF_d1 %>% dplyr::filter(subSystem_CBNR == toString(BD$subSystem[[10]]))
		DF_d2 = DF_d2 %>% dplyr::filter(subSystem_CBNR == toString(BD$subSystem[[10]]))
		DF_d1$disrupted.module <- rep(c(x[[10]][1]),times=nrow(DF_d1))
		DF_d2$disrupted.module <- rep(c(x[[10]][2]),times=nrow(DF_d2))	

		DF_d1 <- subset(DF_d1, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		DF_d2 <- subset(DF_d2, TRUE, c("subSystem_CB", "disrupted.module", "del_rxnID_KO"))
		colnames(DF_d1) = c("subSystem", "disrupted.module", "rxnList")
		colnames(DF_d2) = c("subSystem", "disrupted.module", "rxnList")
		BD_NR_S4 = RBIND(list(DF_d1, DF_d2))

		write.table(BD_NR_S4, "PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S4.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
						
###############################################
# Backtracking - Minerva Visualization
###############################################

	# Make Minerva-friendly file
		
		# Final-subSystem-wise (n=14)
			
			# duplicate BD files
			DF_1 = BD_S1
			DF_2 = BD_S2
			DF_3 = BD_S3
			DF_4 = BD_R_S1
			DF_5 = BD_R_S2
			DF_6 = BD_R_S3
			DF_7 = BD_NR_S1
			DF_8 = BD_NR_S2
			DF_9 = BD_NR_S3
			DF_10 = BD_NR_S4
			
			# add 'r_' as prefix to Rxns
			DF_1$rxnList = paste0('r_', DF_1$rxnList)
			DF_2$rxnList = paste0('r_', DF_2$rxnList)
			DF_3$rxnList = paste0('r_', DF_3$rxnList)
			DF_4$rxnList = paste0('r_', DF_4$rxnList)
			DF_5$rxnList = paste0('r_', DF_5$rxnList)
			DF_6$rxnList = paste0('r_', DF_6$rxnList)
			DF_7$rxnList = paste0('r_', DF_7$rxnList)
			DF_8$rxnList = paste0('r_', DF_8$rxnList)
			DF_9$rxnList = paste0('r_', DF_9$rxnList)
			DF_10$rxnList = paste0('r_', DF_10$rxnList)
						
			# rm duplicates
			DF_1 <- DF_1[!duplicated(DF_1[,c('rxnList')]),]
			DF_2 <- DF_2[!duplicated(DF_2[,c('rxnList')]),]
			DF_3 <- DF_3[!duplicated(DF_3[,c('rxnList')]),]
			DF_4 <- DF_4[!duplicated(DF_4[,c('rxnList')]),]
			DF_5 <- DF_5[!duplicated(DF_5[,c('rxnList')]),]
			DF_6 <- DF_6[!duplicated(DF_6[,c('rxnList')]),]
			DF_7 <- DF_7[!duplicated(DF_7[,c('rxnList')]),]
			DF_8 <- DF_8[!duplicated(DF_8[,c('rxnList')]),]
			DF_9 <- DF_9[!duplicated(DF_9[,c('rxnList')]),]
			DF_10 <- DF_10[!duplicated(DF_10[,c('rxnList')]),]
			
			# add 'lineWidth'
			DF_1$lineWidth <- rep(c("3"),times=nrow(DF_1))
			DF_2$lineWidth <- rep(c("3"),times=nrow(DF_2))
			DF_3$lineWidth <- rep(c("3"),times=nrow(DF_3))
			DF_4$lineWidth <- rep(c("3"),times=nrow(DF_4))
			DF_5$lineWidth <- rep(c("3"),times=nrow(DF_5))
			DF_6$lineWidth <- rep(c("3"),times=nrow(DF_6))
			DF_7$lineWidth <- rep(c("3"),times=nrow(DF_7))
			DF_8$lineWidth <- rep(c("3"),times=nrow(DF_8))
			DF_9$lineWidth <- rep(c("3"),times=nrow(DF_9))
			DF_10$lineWidth <- rep(c("3"),times=nrow(DF_10))
								
			# add 'color'
			DF_1$color <- rep(c("#0000FF"),times=nrow(DF_1))	#blue
			DF_2$color <- rep(c("#0000FF"),times=nrow(DF_2))	#blue
			DF_3$color <- rep(c("#0000FF"),times=nrow(DF_3))	#blue
			DF_4$color <- rep(c("#00FF00"),times=nrow(DF_4))	#green
			DF_5$color <- rep(c("#00FF00"),times=nrow(DF_5))	#green
			DF_6$color <- rep(c("#00FF00"),times=nrow(DF_6))	#green
			DF_7$color <- rep(c("#FF0000"),times=nrow(DF_7))	#red
			DF_8$color <- rep(c("#FF0000"),times=nrow(DF_8))	#red
			DF_9$color <- rep(c("#FF0000"),times=nrow(DF_9))	#red
			DF_10$color <- rep(c("#FF0000"),times=nrow(DF_10))	#red
			
			# slice data
			DF_1 <- subset(DF_1, TRUE, c("rxnList", "lineWidth", "color"))
			DF_2 <- subset(DF_2, TRUE, c("rxnList", "lineWidth", "color"))
			DF_3 <- subset(DF_3, TRUE, c("rxnList", "lineWidth", "color"))
			DF_4 <- subset(DF_4, TRUE, c("rxnList", "lineWidth", "color"))
			DF_5 <- subset(DF_5, TRUE, c("rxnList", "lineWidth", "color"))
			DF_6 <- subset(DF_6, TRUE, c("rxnList", "lineWidth", "color"))
			DF_7 <- subset(DF_7, TRUE, c("rxnList", "lineWidth", "color"))
			DF_8 <- subset(DF_8, TRUE, c("rxnList", "lineWidth", "color"))
			DF_9 <- subset(DF_9, TRUE, c("rxnList", "lineWidth", "color"))
			DF_10 <- subset(DF_10, TRUE, c("rxnList", "lineWidth", "color"))
			
			# rename colnames
			colnames(DF_1) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_2) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_3) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_4) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_5) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_6) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_7) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_8) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_9) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF_10) = c("reactionIdentifier", "lineWidth", "color")
			
			dim(DF_1)
			dim(DF_2)
			dim(DF_3)
			dim(DF_4)
			dim(DF_5)
			dim(DF_6)
			dim(DF_7)
			dim(DF_8)
			dim(DF_9)
			dim(DF_10)
									
			# write table
			write.table(DF_1, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_1_BD_S1_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			write.table(DF_2, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_2_BD_S2_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			write.table(DF_3, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_3_BD_S3_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			write.table(DF_4, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_4_BD_R_S1_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			write.table(DF_5, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_5_BD_R_S2_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			write.table(DF_6, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_6_BD_R_S3_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)		
			write.table(DF_7, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_7_BD_NR_S1_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)		
			write.table(DF_8, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_8_BD_NR_S2_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)		
			write.table(DF_9, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_9_BD_NR_S3_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)		
			write.table(DF_10, "PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_10_BD_NR_S4_Minerva.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)		
		
		# Final-phenotype-wise (n=3)
		
			# Merge DFs for Minerva
			DF1 = RBIND(list(BD_S1, BD_S2, BD_S3))
			DF2 = RBIND(list(BD_R_S1, BD_R_S2, BD_R_S3))
			DF3 = RBIND(list(BD_NR_S1, BD_NR_S2, BD_NR_S3, BD_NR_S4))
			
			# add 'r_' as prefix to Rxns
			DF1$rxnList = paste0('r_', DF1$rxnList)
			DF2$rxnList = paste0('r_', DF2$rxnList)
			DF3$rxnList = paste0('r_', DF3$rxnList)
			
			# rm duplicates
			DF1 <- DF1[!duplicated(DF1[,c('rxnList')]),]
			DF2 <- DF2[!duplicated(DF2[,c('rxnList')]),]
			DF3 <- DF3[!duplicated(DF3[,c('rxnList')]),]
			
			# add 'lineWidth'
			DF1$lineWidth <- rep(c("3"),times=nrow(DF1))
			DF2$lineWidth <- rep(c("3"),times=nrow(DF2))
			DF3$lineWidth <- rep(c("3"),times=nrow(DF3))
					
			# add 'color'
			DF1$color <- rep(c("#0000FF"),times=nrow(DF1))	#blue
			DF2$color <- rep(c("#00FF00"),times=nrow(DF2))	#green
			DF3$color <- rep(c("#FF0000"),times=nrow(DF3))	#red
			
			# slice data
			DF1 <- subset(DF1, TRUE, c("rxnList", "lineWidth", "color"))
			DF2 <- subset(DF2, TRUE, c("rxnList", "lineWidth", "color"))
			DF3 <- subset(DF3, TRUE, c("rxnList", "lineWidth", "color"))
			
			# rename colnames
			colnames(DF1) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF2) = c("reactionIdentifier", "lineWidth", "color")
			colnames(DF3) = c("reactionIdentifier", "lineWidth", "color")
			
#~ 			dim(DF1)
#~ 			dim(DF2)
#~ 			dim(DF3)
			
			# write table
			write.table(DF1, "PlotResults/BD_Rxns.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			write.table(DF2, "PlotResults/BD_R_Rxns.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
			write.table(DF3, "PlotResults/BD_NR_Rxns.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
		
		
		




