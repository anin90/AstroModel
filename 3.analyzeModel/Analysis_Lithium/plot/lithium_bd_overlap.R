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
library(pheatmap)
library(magrittr)
library(venn)

pdf("lithium_bd_hyper.pdf", height=4, width=4)

###########
# Load data
###########
iAstro_Primary_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_Primary_TP.csv", header = T, sep = "\t")
iAstro_iPS_Ctrl_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_Ctrl_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_R_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_R_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_NR_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_NR_TP.csv", header = T, sep = "\t")

bd_lumped_tbl_all <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_TP.csv", header = T, sep = "\t")
bd_r_tbl_all <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_R_TP.csv", header = T, sep = "\t")		
bd_nr_tbl_all <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/FSR_iAstro_iPS_Ctrl_TP_vs_iAstro_iPS_BD_NR_TP.csv", header = T, sep = "\t")

bd_212_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_212_tbl.csv", header = T, sep = "\t")
bd_r_92_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_r_92_tbl.csv", header = T, sep = "\t")		
bd_nr_670_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_nr_670_tbl.csv", header = T, sep = "\t")

bd_lumped_rxns_fdr <- read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_lumped_rxns_fdr.txt", header = T, sep = "\t")
bd_r_rxns_fdr <- read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_r_rxns_fdr.txt", header = T, sep = "\t")		
bd_nr_rxns_fdr <- read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_tbl_significant/bd_nr_rxns_fdr.txt", header = T, sep = "\t")

#~ 	venn(x, ilab=TRUE, zcolor = "style")

###############################################################################################
# https://stats.stackexchange.com/questions/95523/significance-of-overlap-between-multiple-lists
###############################################################################################
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

########################
# Load Primary_TP_Lithium
########################

#~ # Fig.1: (Primary_TP_Lithium vs BD_Disrupted)
#~ 	Akkouh_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_Akkouh_Input.csv",header = T, sep = "\t")
#~ 	GSE66276_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_GSE66276_Input.csv",header = T, sep = "\t")
#~ 	GSE132397_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_GSE132397_Input.csv",header = T, sep = "\t")
#~ 	Akkouh = Akkouh_Input$m_in #InputRxnsOnly
#~ 	GSE66276 = GSE66276_Input$m_in #InputRxnsOnly
#~ 	GSE132397 = GSE132397_Input$m_in #InputRxnsOnly
	
#~ 	#1b (plot_hypermat) - Li_GSE132397_removed
#~ 	x <- list(Li_Akkouh = Akkouh, Li_GSE66276 = GSE66276,	
#~ 	BD_Lumped = bd_212_tbl$rxnList, BD_R = bd_r_92_tbl$rxnList, BD_NR = bd_nr_670_tbl$rxnList)
	
#~ 	y = list(iAstro_Primary_TP$Var1, iAstro_iPS_BD_TP$Var1, iAstro_iPS_BD_R_TP$Var1, iAstro_iPS_BD_NR_TP$Var1);
#~ 	yInt = Reduce(union,y); background = length(yInt);
#~ 	M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
#~ 	M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
#~     rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
#~ 	pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
#~ 	main = 'Li+.rxns vs BD.rxns.all', color=colorRampPalette(c("green", "yellow", "orangered"))(50))
	
#~ 	M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
#~ 	#M.int

#~ # Fig.2: (iPS_Ctrl_TP_Lithium vs BD_Disrupted)
#~ 	Akkouh_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/iPS_Ctrl_TP_Akkouh_Input.csv",header = T, sep = "\t")
#~ 	GSE66276_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/iPS_Ctrl_TP_GSE66276_Input.csv",header = T, sep = "\t")
#~ 	GSE132397_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/iPS_Ctrl_TP_GSE132397_Input.csv",header = T, sep = "\t")
#~ 	Akkouh = Akkouh_Input$m_in #InputRxnsOnly
#~ 	GSE66276 = GSE66276_Input$m_in #InputRxnsOnly
#~ 	GSE132397 = GSE132397_Input$m_in #InputRxnsOnly

#~ 	#2b (plot_hypermat) - Li_GSE132397_removed
#~ 	x <- list(Li_Akkouh = Akkouh, Li_GSE66276 = GSE66276,	
#~ 	BD_Lumped = bd_212_tbl$rxnList, BD_R = bd_r_92_tbl$rxnList, BD_NR = bd_nr_670_tbl$rxnList)
	
#~ 	y = list(iAstro_iPS_Ctrl_TP$Var1, iAstro_iPS_BD_TP$Var1, iAstro_iPS_BD_R_TP$Var1, iAstro_iPS_BD_NR_TP$Var1);
#~ 	yInt = Reduce(union,y); background = length(yInt);
#~ 	M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
#~ 	M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
#~     rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
#~ 	pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
#~ 	main = 'Li+ (iPS_Ctrl) vs BD', color=colorRampPalette(c("green", "yellow", "orangered"))(50))
	
#~ 	M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
#~ 	#M.int

#~ # Fig.3: (iPS_BD_TP_Lithium vs BD_Disrupted)
#~ 	Akkouh_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/iPS_BD_TP_Akkouh_Input.csv",header = T, sep = "\t")
#~ 	GSE66276_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/iPS_BD_TP_GSE66276_Input.csv",header = T, sep = "\t")
#~ 	GSE132397_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/iPS_BD_TP_GSE132397_Input.csv",header = T, sep = "\t")
#~ 	Akkouh = Akkouh_Input$m_in #InputRxnsOnly
#~ 	GSE66276 = GSE66276_Input$m_in #InputRxnsOnly
#~ 	GSE132397 = GSE132397_Input$m_in #InputRxnsOnly
	
#~ 	#3b (plot_hypermat) - Li_GSE132397_removed
#~ 	x <- list(Li_Akkouh = Akkouh, Li_GSE66276 = GSE66276,	
#~ 	BD_Lumped = bd_212_tbl$rxnList, BD_R = bd_r_92_tbl$rxnList, BD_NR = bd_nr_670_tbl$rxnList)
	
#~ 	y = list(iAstro_iPS_BD_TP$Var1, iAstro_iPS_BD_R_TP$Var1, iAstro_iPS_BD_NR_TP$Var1);
#~ 	yInt = Reduce(union,y); background = length(yInt);
#~ 	M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
#~ 	M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
#~     rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
#~ 	pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
#~ 	main = 'Li+ (iPS_BD) vs BD', color=colorRampPalette(c("green", "yellow", "orangered"))(50))

#~ 	M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
#~ 	#M.int
	
#~ ################
#~ # rnxs slice
#~ ################
#~ 	li_bd_Akkouh = intersect(bd_212_tbl$rxnList, Akkouh);	li_bd_r_Akkouh = intersect(bd_r_92_tbl$rxnList, Akkouh);
#~ 	li_bd_nr_Akkouh = intersect(bd_nr_670_tbl$rxnList, Akkouh);	li_bd_GSE66276 = intersect(bd_212_tbl$rxnList, GSE66276);
#~ 	li_bd_r_GSE66276 = intersect(bd_r_92_tbl$rxnList, GSE66276);	li_bd_nr_GSE66276 = intersect(bd_nr_670_tbl$rxnList, GSE66276);
#~ 	x = c(length(li_bd_Akkouh), length(li_bd_r_Akkouh), length(li_bd_nr_Akkouh), length(li_bd_GSE66276), length(li_bd_r_GSE66276), 
#~ 			length(li_bd_nr_GSE66276)); x;


########################
# Load Primary_TP_Lithium
########################

#~ # Fig.5: (Primary_TP_Lithium vs BD_Disrupted_all)
#~ 	Akkouh_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_Akkouh_Input.csv",header = T, sep = "\t")
#~ 	GSE66276_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_GSE66276_Input.csv",header = T, sep = "\t")
#~ 	GSE132397_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_GSE132397_Input.csv",header = T, sep = "\t")
#~ 	Akkouh = Akkouh_Input$m_in #InputRxnsOnly
#~ 	GSE66276 = GSE66276_Input$m_in #InputRxnsOnly
#~ 	GSE132397 = GSE132397_Input$m_in #InputRxnsOnly
	
#~ 	#5a (plot_hypermat) - Li_GSE132397_removed
#~ 	x <- list(Li_Akkouh = Akkouh, Li_GSE66276 = GSE66276,	
#~ 	BD_Lumped = bd_lumped_tbl_all$m, BD_R = bd_r_tbl_all$m, BD_NR = bd_nr_tbl_all$m)
	
#~ 	y = list(iAstro_Primary_TP$Var1, iAstro_iPS_BD_TP$Var1, iAstro_iPS_BD_R_TP$Var1, iAstro_iPS_BD_NR_TP$Var1);
#~ 	yInt = Reduce(intersect,y); background = length(yInt);
#~ 	M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
#~ 	M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
#~     rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
#~ 	pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
#~ 	main = 'Li+.rxns vs BD.rxns.all', color=colorRampPalette(c("green", "yellow", "orangered"))(50))
	
#~ 	M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
#~ 	M.int

# Fig.4: (Primary_TP_Lithium vs BD_Disrupted_fdr)
	Akkouh_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_Akkouh_Input.csv",header = T, sep = "\t")
	GSE66276_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_GSE66276_Input.csv",header = T, sep = "\t")
	GSE132397_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/Primary_TP_GSE132397_Input.csv",header = T, sep = "\t")
	Akkouh = Akkouh_Input$m_in #InputRxnsOnly
	GSE66276 = GSE66276_Input$m_in #InputRxnsOnly
	GSE132397 = GSE132397_Input$m_in #InputRxnsOnly
	
	#4a (plot_hypermat) - Li_GSE132397_removed
	x <- list(Li_Akkouh = Akkouh, Li_GSE66276 = GSE66276,	
	BD_Lumped = bd_lumped_rxns_fdr$rxnList, BD_R = bd_r_rxns_fdr$rxnList, BD_NR = bd_nr_rxns_fdr$rxnList)
	
	y = list(iAstro_Primary_TP$Var1, iAstro_iPS_BD_TP$Var1, iAstro_iPS_BD_R_TP$Var1, iAstro_iPS_BD_NR_TP$Var1);
	yInt = Reduce(intersect,y); background = length(yInt);
	M <- hyper_matrix(x, background);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
	M.adj = M %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(M));
    rownames(M.adj) = rownames(M); colnames(M.adj) = colnames(M);
	pheatmap(M.adj, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, number_format = "%.1e",
	main = 'Li+.rxns vs BD.rxns.fdr', color=colorRampPalette(c("green", "yellow", "orangered"))(50))
	
	M.int <- hyper_matrix(x, background);	M.int[lower.tri(M.int)] <- NA;	
	M.int

	#rnxs slice
	li_bd_fdr_Akkouh = intersect(bd_lumped_rxns_fdr$rxnList, Akkouh);
	li_bd_r_fdr_Akkouh = intersect(bd_r_rxns_fdr$rxnList, Akkouh);
	li_bd_nr_fdr_Akkouh = intersect(bd_nr_rxns_fdr$rxnList, Akkouh);	
	
	li_bd_fdr_GSE66276 = intersect(bd_lumped_rxns_fdr$rxnList, GSE66276);
	li_bd_r_fdr_GSE66276 = intersect(bd_r_rxns_fdr$rxnList, GSE66276);
	li_bd_nr_fdr_GSE66276 = intersect(bd_nr_rxns_fdr$rxnList, GSE66276);
	
	x = c(length(li_bd_fdr_Akkouh), length(li_bd_r_fdr_Akkouh), length(li_bd_nr_fdr_Akkouh),
			length(li_bd_fdr_GSE66276), length(li_bd_r_fdr_GSE66276), length(li_bd_nr_fdr_GSE66276)); x;

	#write fdr.significant subset to csv		
	#Li_Akkouh_vs_BD_Lumped_FDR	
	Li_Akkouh_vs_BD_Lumped_FDR_rxns = bd_lumped_rxns_fdr[bd_lumped_rxns_fdr$rxnList %in% li_bd_fdr_Akkouh,]
	write.table(Li_Akkouh_vs_BD_Lumped_FDR_rxns, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_Akkouh_vs_BD_Lumped_FDR_rxns.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	Li_Akkouh_vs_BD_Lumped_FDR_subSystem = Li_Akkouh_vs_BD_Lumped_FDR_rxns %>% count(subSystem)
	write.table(Li_Akkouh_vs_BD_Lumped_FDR_subSystem, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_Akkouh_vs_BD_Lumped_FDR_subSystem.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	head(Li_Akkouh_vs_BD_Lumped_FDR_subSystem)
	#Li_Akkouh_vs_BD_R_FDR	
	Li_Akkouh_vs_BD_R_FDR_rxns = bd_r_rxns_fdr[bd_r_rxns_fdr$rxnList %in% li_bd_r_fdr_Akkouh,]
	write.table(Li_Akkouh_vs_BD_R_FDR_rxns, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_Akkouh_vs_BD_R_FDR_rxns.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	Li_Akkouh_vs_BD_R_FDR_subSystem = Li_Akkouh_vs_BD_R_FDR_rxns %>% count(subSystem)
	write.table(Li_Akkouh_vs_BD_R_FDR_subSystem, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_Akkouh_vs_BD_R_FDR_subSystem.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	head(Li_Akkouh_vs_BD_R_FDR_subSystem)
	#Li_Akkouh_vs_BD_NR_FDR	
	Li_Akkouh_vs_BD_NR_FDR_rxns = bd_nr_rxns_fdr[bd_nr_rxns_fdr$rxnList %in% li_bd_nr_fdr_Akkouh,]
	write.table(Li_Akkouh_vs_BD_NR_FDR_rxns, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_Akkouh_vs_BD_NR_FDR_rxns.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	Li_Akkouh_vs_BD_NR_FDR_subSystem = Li_Akkouh_vs_BD_NR_FDR_rxns %>% count(subSystem)
	write.table(Li_Akkouh_vs_BD_NR_FDR_subSystem, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_Akkouh_vs_BD_NR_FDR_subSystem.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	head(Li_Akkouh_vs_BD_NR_FDR_subSystem)
	
	#Li_GSE66276_vs_BD_Lumped_FDR	
	Li_GSE66276_vs_BD_Lumped_FDR_rxns = bd_lumped_rxns_fdr[bd_lumped_rxns_fdr$rxnList %in% li_bd_fdr_GSE66276,]
	write.table(Li_GSE66276_vs_BD_Lumped_FDR_rxns, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_GSE66276_vs_BD_Lumped_FDR_rxns.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	Li_GSE66276_vs_BD_Lumped_FDR_subSystem = Li_GSE66276_vs_BD_Lumped_FDR_rxns %>% count(subSystem)
	write.table(Li_GSE66276_vs_BD_Lumped_FDR_subSystem, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_GSE66276_vs_BD_Lumped_FDR_subSystem.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	head(Li_GSE66276_vs_BD_Lumped_FDR_subSystem)
	#Li_GSE66276_vs_BD_R_FDR	
	Li_GSE66276_vs_BD_R_FDR_rxns = bd_r_rxns_fdr[bd_r_rxns_fdr$rxnList %in% li_bd_r_fdr_GSE66276,]
	write.table(Li_GSE66276_vs_BD_R_FDR_rxns, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_GSE66276_vs_BD_R_FDR_rxns.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	Li_GSE66276_vs_BD_R_FDR_subSystem = Li_GSE66276_vs_BD_R_FDR_rxns %>% count(subSystem)
	write.table(Li_GSE66276_vs_BD_R_FDR_subSystem, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_GSE66276_vs_BD_R_FDR_subSystem.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	head(Li_GSE66276_vs_BD_R_FDR_subSystem)
	#Li_GSE66276_vs_BD_NR_FDR	
	Li_GSE66276_vs_BD_NR_FDR_rxns = bd_nr_rxns_fdr[bd_nr_rxns_fdr$rxnList %in% li_bd_nr_fdr_GSE66276,]
	write.table(Li_GSE66276_vs_BD_NR_FDR_rxns, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_GSE66276_vs_BD_NR_FDR_rxns.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	Li_GSE66276_vs_BD_NR_FDR_subSystem = Li_GSE66276_vs_BD_NR_FDR_rxns %>% count(subSystem)
	write.table(Li_GSE66276_vs_BD_NR_FDR_subSystem, "/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/plot/Li_vs_BD_FDR_Overlap/Li_GSE66276_vs_BD_NR_FDR_subSystem.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	head(Li_GSE66276_vs_BD_NR_FDR_subSystem)
	

