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

pdf("lithium_bd_venn.pdf")

###########
# Load data
###########
bd_212_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_212_tbl.csv", header = T, sep = "\t")
bd_r_92_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_r_92_tbl.csv", header = T, sep = "\t")		
bd_nr_670_tbl <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_BD/plot/bd_nr_670_tbl.csv", header = T, sep = "\t")
iAstro_Primary_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_Primary_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_R_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_R_TP.csv", header = T, sep = "\t")
iAstro_iPS_BD_NR_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Rxns_iAstro_iPS_BD_NR_TP.csv", header = T, sep = "\t")

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
  # round to 2 digit.
  return(round(combination,2))
}

########################
# Load Primary_TP_Akkouh
########################
Untreated_vs_Lithium_Output <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/iAstro_Primary_TP_Untreated_vs_Lithium_Output.csv",header = T, sep = "\t")
Untreated_vs_Lithium_Input <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_Lithium/iAstro_Primary_TP_Untreated_vs_Lithium_Input.csv",header = T, sep = "\t")
bind = bind_rows(Untreated_vs_Lithium_Output,Untreated_vs_Lithium_Input)
Untreated_vs_Lithium_all = coalesce(bind$m_out, bind$m_in)

# Fig.1: (Primary_TP_Akkouh vs BD_Disrupted)

#1a (venn rxns)
x <- list(
  Akkouh = Untreated_vs_Lithium_all,	
  BD_Lumped = bd_212_tbl$rxnList, 
  BD_R = bd_r_92_tbl$rxnList, 
  BD_NR = bd_nr_670_tbl$rxnList)
p = ggvenn(x, fill_color = c("grey","blue", "green", "red")); p
#svg(filename="lithium_bd_venn.svg"); p; dev.off()

#1b (venn hypergeometric test)
background <- sum(c(nrow(Untreated_vs_Lithium_all),nrow(bd_212_tbl), nrow(bd_r_92_tbl), nrow(bd_nr_670_tbl)))				
M <- hyper_matrix(x, 6500);	M[upper.tri(M)] <- NA;	diag(M)<-NA;
pheatmap(M, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE, main = 'Significance of pairwise overlap', 
	fontsize_number = 20, fontsize = 15)
#p = pheatmap(M, cluster_rows=F, cluster_cols=F, na_col="white", display_numbers = TRUE)	
#svg(filename="lithium_bd_hypermat.svg", width=3, height=3); p; dev.off()


