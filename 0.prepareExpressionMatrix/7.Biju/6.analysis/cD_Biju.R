#############
start.time <- Sys.time()
#############
#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
pacman::p_load(affy, pheatmap, RColorBrewer, dplyr, tidyverse, annotate, rat2302.db, mouse4302.db, homologene, 
readxl, oligo, limma, mogene20sttranscriptcluster.db, qvalue, GEOquery, tidyr, tibble, splitstackshape, gplots, 
ggplot2, ggfortify, reshape2, factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, 
mgsub, gsubfn, readxl, openxlsx, UpSetR, qvalue, GEOquery, TeachingDemos, sm, org.Hs.eg.db, data.table) 

setwd("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/0.prepareExpressionMatrix/7.Biju/6.analysis")

#############
# Load data
#############

	fpkm_redone <- read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/0.prepareExpressionMatrix/7.Biju/6.analysis/genes.fpkm_tracking",sep="\t",header=TRUE)
	fpkm_redone_subset = fpkm_redone[,!grepl("*conf",names(fpkm_redone))] #drop unnessaary fields
	fpkm_redone_subset = fpkm_redone[,!grepl("*status",names(fpkm_redone))] #drop unnessaary fields

	# slice.data
	fpkm_redone_subset = subset(fpkm_redone_subset, TRUE, c("gene_short_name","q1_FPKM","q2_FPKM","q3_FPKM","q4_FPKM","q5_FPKM","q6_FPKM","q7_FPKM","q8_FPKM",
	"q9_FPKM","q10_FPKM","q11_FPKM","q12_FPKM","q13_FPKM","q14_FPKM","q15_FPKM","q16_FPKM","q17_FPKM","q18_FPKM","q19_FPKM","q20_FPKM","q21_FPKM","q22_FPKM","q23_FPKM","q24_FPKM"))

	# rename.cols.level.1
	colnames(fpkm_redone_subset) = c("gene_redone","20007A_S1","20007B_S2","20008A_S3","20008B_S4","20009A_S5","20009B_S6",
	"20010A_S7","20010B_S8","20011A_S9","20011B_S10","20012A_S11","20012B_S12","20013A_S13","20013B_S26","20014A_S15",
	"20014B_S16","20015A_S17","20015B_S18","20016A_S19","20016B_S20","20017A_S21","20017B_S22","20018A_S23","20018B_S24")

	# rename.cols.level.2
	colnames(fpkm_redone_subset) = c("gene_redone","14472.Veh.R1","14470.Veh.R1","14594.Veh.R1","NIH1.Veh.R1","14472.Li.R1","14470.Li.R1",
	"14594.Li.R1","NIH1.Li.R1","14472.Val.R1","14470.Val.R1","14594.Val.R1","NIH1.Val.R1","14472.Veh.R2","14470.Veh.R2","14594.Veh.R2",
	"NIH1.Veh.R2","14472.Li.R2","14470.Li.R2","NIH1.Li.R2","14472.Val.R2","14470.Val.R2","14594.Val.R2","NIH1.Val.R2","14594.Li.R2")

	# rename.cols.level.3
	colnames(fpkm_redone_subset) = c("gene_redone","BD_NR.Veh.R1","BD_R.Veh.R1","Ctrl_Fam.Veh.R1","Ctrl_NIH1.Veh.R1","BD_NR.Li.R1","BD_R.Li.R1",
	"Ctrl_Fam.Li.R1","Ctrl_NIH1.Li.R1","BD_NR.Val.R1","BD_R.Val.R1","Ctrl_Fam.Val.R1","Ctrl_NIH1.Val.R1","BD_NR.Veh.R2","BD_R.Veh.R2","Ctrl_Fam.Veh.R2",
	"Ctrl_NIH1.Veh.R2","BD_NR.Li.R2","BD_R.Li.R2","Ctrl_NIH1.Li.R2","BD_NR.Val.R2","BD_R.Val.R2","Ctrl_Fam.Val.R2","Ctrl_NIH1.Val.R2","Ctrl_Fam.Li.R2")

#############################
# Replace rownames by index
#############################

	tt = cbind(row.names(fpkm_redone_subset), fpkm_redone_subset)
	rownames(tt) = NULL

#############################
# Annotate gene symbols
#############################

	cols <- c("ENTREZID", "SYMBOL")
	gene_list <- as.vector(unlist(tt$gene_redone))
	df = select(org.Hs.eg.db, keys=gene_list, columns=cols, keytype="SYMBOL")

	# remove non-specific gene_symbols (same symbols mapping to multiple entrezids)
	df_clean_probe = df %>% 
			group_by(SYMBOL) %>%  ## group by id column 
			filter(!any(row_number() > 1))

	# remove non-specific entrez-IDs  (multiple probes mapping to the same gene
	df_clean_entrezid = df_clean_probe %>% 
			group_by(ENTREZID) %>%  ## group by id column 
			filter(!any(row_number() > 1))

	# returns TRUE, if there are no duplicates -> before merging
	length(unique(df_clean_entrezid$SYMBOL)) == nrow(df_clean_entrezid)
	length(unique(df_clean_entrezid$ENTREZID)) == nrow(df_clean_entrezid)	

	# save cleaned data frame
	df_clean = df_clean_entrezid

	# merge annotation with original data
	comb = merge(df_clean, tt, by.x = "SYMBOL", by.y = "gene_redone")

	# remove lines that contain NAs across 'ENTREZID' & 'SYMBOL' columns
	fpkm_redone_subset_anno = comb[complete.cases(comb), ]

	# returns TRUE, if there are no duplicates -> after merging
	length(unique(fpkm_redone_subset_anno$SYMBOL)) == nrow(fpkm_redone_subset_anno)
	length(unique(fpkm_redone_subset_anno$ENTREZID)) == nrow(fpkm_redone_subset_anno)

	# change the row name to ENTREZID
	row.names(fpkm_redone_subset_anno) <- fpkm_redone_subset_anno$SYMBOL

	# retain genes expressed (fpkm>=0.1) in atleast 50% (n=12) of all samples (n=24)
	fpkm_redone_subset_anno_filt = fpkm_redone_subset_anno %>% filter(rowSums(.[4:ncol(fpkm_redone_subset_anno)]>=0.1)>=12)
	fpkm_redone_subset_anno_filt = fpkm_redone_subset_anno_filt[ , -which(names(fpkm_redone_subset_anno_filt) %in% c("row.names(fpkm_redone_subset)"))]

	# write final output to table
	write.table(fpkm_redone_subset_anno_filt, "fpkm_redone_subset_anno_filt.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)

	# genes expressed with FPKM >= 0.1 in each sample
	dat = colSums(fpkm_redone_subset_anno_filt[,3:ncol(fpkm_redone_subset_anno_filt)]>=0.1)
	dotchart(dat,xlab="#Genes with FPKM >= 0.1",cex=.7)

	# subset only Veh, Li and Val
		# Veh
		fpkm_redone_subset_anno_filt_Ctrl_NIH1_Veh = subset(fpkm_redone_subset_anno_filt, TRUE, c("SYMBOL", "ENTREZID",
											"Ctrl_NIH1.Veh.R1","Ctrl_NIH1.Veh.R2")) #Ctrl_NIH
		fpkm_redone_subset_anno_filt_Ctrl_Fam_Veh = subset(fpkm_redone_subset_anno_filt, TRUE, c("SYMBOL", "ENTREZID",
											"Ctrl_Fam.Veh.R1","Ctrl_Fam.Veh.R2")) #Ctrl_Fam
		fpkm_redone_subset_anno_filt_BD_R_Fam_Veh = subset(fpkm_redone_subset_anno_filt, TRUE, c("SYMBOL", "ENTREZID",
											"BD_R.Veh.R1","BD_R.Veh.R2")) #BD_R
		fpkm_redone_subset_anno_filt_BD_NR_Fam_Veh = subset(fpkm_redone_subset_anno_filt, TRUE, c("SYMBOL", "ENTREZID",
											"BD_NR.Veh.R1","BD_NR.Veh.R2"))	 #BD_NR							
		fpkm_redone_subset_anno_filt_BD_Fam_Veh = fpkm_redone_subset_anno_filt[,!grepl("*Ctrl|*Li|*Val",
											names(fpkm_redone_subset_anno_filt))] #BD_lumped

		# Li (TBD)
		# Val (TBD)
											
	# write final output_subset to table
	write.table(fpkm_redone_subset_anno_filt_Ctrl_NIH1_Veh, "fpkm_redone_subset_anno_filt_Ctrl_NIH1_Veh.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	write.table(fpkm_redone_subset_anno_filt_Ctrl_Fam_Veh, "fpkm_redone_subset_anno_filt_Ctrl_Fam_Veh.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	write.table(fpkm_redone_subset_anno_filt_BD_R_Fam_Veh, "fpkm_redone_subset_anno_filt_BD_R_Fam_Veh.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	write.table(fpkm_redone_subset_anno_filt_BD_NR_Fam_Veh, "fpkm_redone_subset_anno_filt_BD_NR_Fam_Veh.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
	write.table(fpkm_redone_subset_anno_filt_BD_Fam_Veh, "fpkm_redone_subset_anno_filt_BD_Fam_Veh.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
					
	dim(fpkm_redone_subset)
	dim(comb)				
	dim(fpkm_redone_subset_anno)
	dim(fpkm_redone_subset_anno_filt)

#############
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#############








