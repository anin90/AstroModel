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

		write.table(mat, "PlotResults/Tbl_abs.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/Tbl_abs.csv", header = T, sep = "\t")

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

		write.table(mat, "PlotResults/Tbl_norm_t1.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/Tbl_norm_t1.csv", header = T, sep = "\t")

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

		write.table(mat, "PlotResults/Tbl_norm_t2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/Tbl_norm_t2.csv", header = T, sep = "\t")

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

		write.table(mat, "PlotResults/Tbl_all.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/Tbl_all.csv", header = T, sep = "\t")

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

		write.table(mat, "PlotResults/Tbl_bd.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/Tbl_bd.csv", header = T, sep = "\t")

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
		write.table(mat_keep, "PlotResults/Tbl_bd_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)	
		
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

		write.table(mat, "PlotResults/Tbl_bd_r.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/Tbl_bd_r.csv", header = T, sep = "\t")

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
		write.table(mat_keep, "PlotResults/Tbl_bd_r_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)	

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

		write.table(mat, "PlotResults/Tbl_bd_nr.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/Tbl_bd_nr.csv", header = T, sep = "\t")

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
		write.table(mat_keep, "PlotResults/Tbl_bd_nr_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
		
		
###########################
# subSystem - Module matrix
###########################

	# tbl_filt
		BD_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/Tbl_bd_filt.csv", header = T, sep = "\t")
		BD_R_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/Tbl_bd_r_filt.csv", header = T, sep = "\t")
		BD_NR_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/PlotResults/Tbl_bd_nr_filt.csv", header = T, sep = "\t")
		
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
		
		write.table(mat, "PlotResults/Tbl_mat_filt.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)
		
		mat <- read.csv("PlotResults/Tbl_mat_filt.csv", header = T, sep = "\t")	
		
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
		write.table(BD_S1, "PlotResults/BD_S1.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
		

