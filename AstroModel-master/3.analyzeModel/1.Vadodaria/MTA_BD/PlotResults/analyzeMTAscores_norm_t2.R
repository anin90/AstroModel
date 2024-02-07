#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
pacman::p_load(pheatmap, RColorBrewer, dplyr, tidyverse, tidyr, tibble, splitstackshape, gplots, ggplot2, ggfortify, reshape2, 
factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, mgsub, gsubfn, readxl, openxlsx, UpSetR) 

pdf("mta_tbl_prctile_top_norm_t2/analyzeMTAscores_norm_t2.pdf")

###################
#~ Load covid_mta
###################  

# Specifying the path name
#~ 	path <- "/home/anirudh/Downloads/papers_to_read/msb202110260-sup-0007-tableev6.xlsx"
#~ 	sheets <- openxlsx::getSheetNames(path)
#~ 	data_frame <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=path)
#~ 	names(data_frame) <- sheets

# Load sheets
#~ 	NHBE = data_frame$"NHBE"
#~ 	Calu_3 = data_frame$"Calu-3"
#~ 	A549 = data_frame$"A549"
#~ 	Vero = data_frame$"Vero"
#~ 	T293 = data_frame$"293T"
#~ 	Caco_2 = data_frame$"Caco-2"
#~ 	Swab.Butler = data_frame$"Swab.Butler"
#~ 	Swab.Lieberman = data_frame$"Swab.Lieberman"
#~ 	BALF = data_frame$"BALF"
#~ 	SC.Liao = data_frame$"SC.Liao"
#~ 	SC.Chua.Basal = data_frame$"SC.Chua.Basal"
#~ 	SC.Chua.Ciliated = data_frame$"SC.Chua.Ciliated"

# Histogram mta_score
#~ 	par(mfrow=c(4,3))
#~ 	hist(NHBE$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="NHBE")
#~ 	hist(Calu_3$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Calu-3")
#~ 	hist(A549$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="A549")
#~ 	hist(Vero$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Vero")
#~ 	hist(T293$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="293T")
#~ 	hist(Caco_2$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Caco-2")
#~ 	hist(Swab.Butler$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Swab.Butler")
#~ 	hist(Swab.Lieberman$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Swab.Lieberman")
#~ 	hist(BALF$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="BALF")
#~ 	hist(SC.Liao$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="SC.Liao")
#~ 	hist(SC.Chua.Basal$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="SC.Chua.Basal")
#~ 	hist(SC.Chua.Ciliated$score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="SC.Chua.Ciliated")

#######################
#~ Load vadodaria_mta
#######################

# Load data
	Ctrl_to_Primary = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_Ctrl_to_Primary.csv", header = T, sep = "\t")
	Primary_to_Ctrl = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_Primary_to_Ctrl.csv", header = T, sep = "\t")
	Ctrl_to_BD = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_Ctrl_to_BD.csv", header = T, sep = "\t")
	BD_to_Ctrl = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_BD_to_Ctrl.csv", header = T, sep = "\t")
	Ctrl_to_BD_R = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_Ctrl_to_BD_R.csv", header = T, sep = "\t")
	BD_R_to_Ctrl = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_BD_R_to_Ctrl.csv", header = T, sep = "\t")
	Ctrl_to_BD_NR = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_Ctrl_to_BD_NR.csv", header = T, sep = "\t")
	BD_NR_to_Ctrl = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_BD_NR_to_Ctrl.csv", header = T, sep = "\t")
	R_to_NR = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_R_to_NR.csv", header = T, sep = "\t")
	NR_to_R = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_norm_t2/Vadodaria_NR_to_R.csv", header = T, sep = "\t")

# Rename data
	CP = Ctrl_to_Primary
	PC = Primary_to_Ctrl
	CB = Ctrl_to_BD
	BC = BD_to_Ctrl
	CBR = Ctrl_to_BD_R
	BRC = BD_R_to_Ctrl
	CBNR = Ctrl_to_BD_NR
	BNRC = BD_NR_to_Ctrl
	RNR = R_to_NR
	NRR = NR_to_R

# Histogram mta_score
	par(mfrow=c(4,2))
	hist(CP$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Ctrl_to_Primary")
	hist(PC$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Primary_to_Ctrl")
	hist(CB$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Ctrl_to_BD")
	hist(BC$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="BD_to_Ctrl")
	hist(CBR$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Ctrl_to_BD_R")
	hist(BRC$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="BD_R_to_Ctrl")
	hist(CBNR$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Ctrl_to_BD_NR")
	hist(BNRC$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="BD_NR_to_Ctrl")
	hist(RNR$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="R_to_NR")
	hist(NRR$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="NR_to_R")

# FilterMTAscore, all
	FilterMTAscore <- function(x) {
		mta_score_prctile <- ecdf(x$mta_score)(x$mta_score)
		mta_score_prctile <- mta_score_prctile*100
		x <- cbind(x, mta_score_prctile)
		x <- arrange(x, -mta_score_prctile)
		return(x) 
		}		

	dat = list(CP, PC, CB, BC, CBR, BRC, CBNR, BNRC, RNR, NRR)
	dat_all = lapply(dat, FilterMTAscore)
	
	CP_prctile = dat_all[[1]]
	PC_prctile = dat_all[[2]]
	CB_prctile = dat_all[[3]]
	BC_prctile = dat_all[[4]]
	CBR_prctile = dat_all[[5]]
	BRC_prctile = dat_all[[6]]
	CBNR_prctile = dat_all[[7]]
	BNRC_prctile = dat_all[[8]]
	RNR_prctile = dat_all[[9]]
	NRR_prctile = dat_all[[10]]
	
	write.table(dat_all[[1]], "mta_tbl_prctile_all_norm_t2/Ctrl_to_Primary.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[2]], "mta_tbl_prctile_all_norm_t2/Primary_to_Ctrl.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[3]], "mta_tbl_prctile_all_norm_t2/Ctrl_to_BD.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[4]], "mta_tbl_prctile_all_norm_t2/BD_to_Ctrl.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[5]], "mta_tbl_prctile_all_norm_t2/Ctrl_to_BD_R.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[6]], "mta_tbl_prctile_all_norm_t2/BD_R_to_Ctrl.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[7]], "mta_tbl_prctile_all_norm_t2/Ctrl_to_BD_NR.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[8]], "mta_tbl_prctile_all_norm_t2/BD_NR_to_Ctrl.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[9]], "mta_tbl_prctile_all_norm_t2/R_to_NR.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[10]], "mta_tbl_prctile_all_norm_t2/NR_to_R.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)


# FilterMTAscore, bottom 80% (unchanged)
	FilterMTAscore_unchanged <- function(x) {
		mta_score_prctile <- ecdf(x$mta_score)(x$mta_score)
		mta_score_prctile <- mta_score_prctile*100
		x <- cbind(x, mta_score_prctile)
		x <- arrange(x, -mta_score_prctile)
		x <- filter(x, mta_score_prctile < 80)
		return(x) 
		}		
	dat = list(CP, PC)
	dat_unchanged = lapply(dat, FilterMTAscore_unchanged)

	CP_filt = dat_unchanged[[1]]
	PC_filt = dat_unchanged[[2]]

# FilterMTAscore, top 20% (significant)
	FilterMTAscore_significant <- function(x) {
		mta_score_prctile <- ecdf(x$mta_score)(x$mta_score)
		mta_score_prctile <- mta_score_prctile*100
		x <- cbind(x, mta_score_prctile)
		x <- arrange(x, -mta_score_prctile)
		x <- filter(x, mta_score_prctile >= 80)
		return(x) 
		}	
	dat = list(CB, BC, CBR, BRC, CBNR, BNRC, RNR, NRR)			
	dat_significant = lapply(dat, FilterMTAscore_significant)

	CB_filt = dat_significant[[1]]
	BC_filt = dat_significant[[2]]
	CBR_filt = dat_significant[[3]]
	BRC_filt = dat_significant[[4]]
	CBNR_filt = dat_significant[[5]]
	BNRC_filt = dat_significant[[6]]
	RNR_filt = dat_significant[[7]]
	NRR_filt = dat_significant[[8]]
	
##################
# Load MTA rxns
##################

	CP_rxns = CP_filt$del_rxnID_KO
	PC_rxns = PC_filt$del_rxnID_KO
	CB_rxns = CB_filt$del_rxnID_KO
	BC_rxns = BC_filt$del_rxnID_KO
	CBR_rxns = CBR_filt$del_rxnID_KO
	BRC_rxns = BRC_filt$del_rxnID_KO
	CBNR_rxns = CBNR_filt$del_rxnID_KO
	BNRC_rxns = BNRC_filt$del_rxnID_KO
	RNR_rxns = RNR_filt$del_rxnID_KO
	NRR_rxns = NRR_filt$del_rxnID_KO

################
# rnxs slice
################

	Ctrl_unchanged = intersect(CP_rxns, PC_rxns)
	BD_top = union(CB_rxns, BC_rxns)
	BD_R_top = union(CBR_rxns, BRC_rxns)
	BD_NR_top = union(CBNR_rxns, BNRC_rxns)
	#~ BD_Response_rxns = union(RNR_rxns, RNR_rxns)

###################################
# rnxs to matrix - wrangle
###################################

	mat = lst(Ctrl_unchanged, BD_top, BD_R_top, BD_NR_top) %>% 
				
	  enframe %>% 
	  unnest %>% 
	  count(name, value) %>% 
	  spread(value, n, fill = 0)

	mat = t(mat)

	write.table(mat, "mta_tbl_prctile_top_norm_t2/MTA_Table_BD_norm_t2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

	mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/PlotResults/mta_tbl_prctile_top_norm_t2/MTA_Table_BD_norm_t2.csv",
			header = T, sep = "\t")

##################################
# rnxs to matrix - UpSetR (TP)
##################################

	p = upset(mat,
		order.by = c("degree","freq"), #query.legend = "bottom",
		matrix.color="black", sets = c("Ctrl_unchanged", "BD_top", "BD_R_top", "BD_NR_top"), 
		keep.order = TRUE, 
		queries = list(list(query = intersects, params = list("BD_NR_top", "BD_R_top", "BD_top", "Ctrl_unchanged"), 
		color = "blue", active = T), 
		list(query = intersects, params = list("BD_R_top", "BD_top", "Ctrl_unchanged"), color = "blue", active = T),
		list(query = intersects, params = list("BD_NR_top", "BD_top", "Ctrl_unchanged"), color = "blue", active = T),
		list(query = intersects, params = list("BD_top", "Ctrl_unchanged"), color = "blue", active = T),
		list(query = intersects, params = list("BD_R_top", "Ctrl_unchanged"), color = "green", active = T),
		list(query = intersects, params = list("BD_NR_top", "Ctrl_unchanged"), color = "red", active = T)),
		
		sets.bar.color=c("dark grey", "dark blue", "dark blue", "dark blue"), nintersects = 60)

	p

################
# rnxs slice
################

	#bd_lumped_a
	x = list(Ctrl_unchanged, BD_top, BD_R_top, BD_NR_top)
	bd_lumped_a = Reduce(intersect, x)		
	length(bd_lumped_a)
	#~ write.table(bd_lumped_a, "mta_tbl_prctile_top_norm_t2/bd_lumped_a.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_lumped_b
	x = list(Ctrl_unchanged, BD_top, BD_R_top)
	xInt = Reduce(intersect, x)
	y = BD_NR_top
	bd_lumped_b = setdiff(xInt,y)
	length(bd_lumped_b)
	#~ write.table(bd_lumped_b, "mta_tbl_prctile_top_norm_t2/bd_lumped_b.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_lumped_c
	x = list(Ctrl_unchanged, BD_top, BD_NR_top)
	xInt = Reduce(intersect, x)
	y = BD_R_top
	bd_lumped_c = setdiff(xInt,y)
	length(bd_lumped_c)
	#~ write.table(bd_lumped_c, "mta_tbl_prctile_top_norm_t2/bd_lumped_c.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_lumped_d
	x = list(Ctrl_unchanged, BD_top)
	xInt = Reduce(intersect, x)
	y = list(BD_R_top, BD_NR_top)
	yInt = Reduce(union,y)
	bd_lumped_d = setdiff(xInt,yInt)
	length(bd_lumped_d)
	#~ write.table(bd_lumped_d, "mta_tbl_prctile_top_norm_t2/bd_lumped_d.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_r
	x = list(Ctrl_unchanged, BD_R_top)
	xInt = Reduce(intersect, x)
	y = list(BD_top, BD_NR_top)
	yInt = Reduce(union,y)
	bd_r = setdiff(xInt,yInt)
	length(bd_r)
	#~ write.table(bd_r, "mta_tbl_prctile_top_norm_t2/bd_r.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_nr
	x = list(Ctrl_unchanged, BD_NR_top)
	xInt = Reduce(intersect, x)
	y = list(BD_top, BD_R_top)
	yInt = Reduce(union,y)
	bd_nr = setdiff(xInt,yInt)
	length(bd_nr)
	#~ write.table(bd_nr, "mta_tbl_prctile_top_norm_t2/bd_nr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)


################################
# grep prctile for sliced rxns
################################

	#all_rxns
	df_list = list(CP_prctile, PC_prctile, CB_prctile, BC_prctile, 
					CBR_prctile, BRC_prctile, CBNR_prctile, BNRC_prctile)
					
	df_list = df_list %>% reduce(full_join, by='del_rxnID_KO')
	df_list_prctile = df_list[,grepl("*mta_score_prctile",names(df_list))]
	df_list_subSystem = df_list[,grepl("*subSystem",names(df_list))]
	df_list_MetabolicUnits = df_list[,grepl("*MetabolicUnits",names(df_list))]
	df_list_Localization = df_list[,grepl("*Localization",names(df_list))]
	df_list_subset = cbind(df_list$del_rxnID_KO, df_list_prctile, 
							df_list_subSystem, df_list_MetabolicUnits, df_list_Localization)
	
	colnames(df_list_subset) = c("del_rxnID_KO", "prctile_CP", "prctile_PC", "prctile_CB", "prctile_BC", 
					"prctile_CBR", "prctile_BRC", "prctile_CBNR", "prctile_BNRC", 
					"subSystem_CP", "subSystem_PC", "subSystem_CB", "subSystem_BC", 
					"subSystem_CBR", "subSystem_BRC", "subSystem_CBNR", "subSystem_BNRC",
					"MetabolicUnits_CP", "MetabolicUnits_PC", "MetabolicUnits_CB", "MetabolicUnits_BC", 
					"MetabolicUnits_CBR", "MetabolicUnits_BRC", "MetabolicUnits_CBNR", "MetabolicUnits_BNRC",
					"Localization_CP", "Localization_PC", "Localization_CB", "Localization_BC", 
					"Localization_CBR", "Localization_BRC", "Localization_CBNR", "Localization_BNRC")
					
	prctile_tbl_all = df_list_subset

	#bd_lumped
	x = list(bd_lumped_a, bd_lumped_b, bd_lumped_c, bd_lumped_d)
	bd_lumped = Reduce(union,x)
	bd_lumped = prctile_tbl_all[prctile_tbl_all$del_rxnID_KO %in% bd_lumped ,]
	write.table(bd_lumped, "mta_tbl_prctile_top_norm_t2/bd_lumped.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	
	#bd_responder
	x = list(bd_r, bd_lumped_b)
	bd_responder = Reduce(union,x)
	bd_responder = prctile_tbl_all[prctile_tbl_all$del_rxnID_KO %in% bd_responder ,]
	write.table(bd_responder, "mta_tbl_prctile_top_norm_t2/bd_responder.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	
	#bd_nonresponder
	x = list(bd_nr, bd_lumped_c)
	bd_nonresponder = Reduce(union,x)
	bd_nonresponder = prctile_tbl_all[prctile_tbl_all$del_rxnID_KO %in% bd_nonresponder ,]
	write.table(bd_nonresponder, "mta_tbl_prctile_top_norm_t2/bd_nonresponder.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)

###########
# Plot data
###########
		
# Fig.1: venn rxns
	x <- list(
		BD_Lumped = bd_lumped$del_rxnID_KO, 
		BD_Responder = bd_responder$del_rxnID_KO, 
		BD_NonResponder = bd_nonresponder$del_rxnID_KO)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlap in rxns disrupted between models")
		
# Fig.2: venn subsystems
	x <- list(
		BD_Lumped = bd_lumped$subSystem_BC, 
		BD_Responder = bd_responder$subSystem_BRC, 
		BD_NonResponder = bd_nonresponder$subSystem_BNRC)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlap in subSystems disrupted between models")		
	
# Fig.3: Number of rxns per subSystem (bd_lumped)
	par(mfrow=c(1,1))
	iPS_BD_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_TP_norm_t2.csv", header = T, sep = "\t")
	iPS_BD_TP_subSystem = iPS_BD_TP %>% count(SubSystem)
	bd_lumped_subSystem = bd_lumped %>% count(subSystem_BC)
	comb = merge(bd_lumped_subSystem, iPS_BD_TP_subSystem, by.x = "subSystem_BC", by.y = "SubSystem")
	comb = comb[order(comb$n.x),]
	#3a
	dotchart(comb$n.x, labels = comb$subSystem_BC,
			cex = 0.5, bg = "blue", xlab = "# of rxns", main = "# of mta hits (bd_lumped)")
	#3b
	dotchart(comb$n.y, labels = comb$subSystem_BC,
			cex = 0.5, bg = "grey", xlab = "# of rxns", main = "# of mta hits (bd_lumped) vs all rxns")	
	points(comb$n.x, 1:nrow(comb), col = "blue", cex = 0.7)	
	#3c
	p.val = phyper(q = comb$n.x-1, m = comb$n.y, n = nrow(iPS_BD_TP)-comb$n.y, k = sum(comb$n.x), lower.tail = FALSE, log.p = FALSE)
	p.val.fdr = p.adjust(p.val, method = "fdr")
	combStat = cbind(comb, p.val, p.val.fdr)
	combStat = combStat[order(-combStat$p.val.fdr),]
	combStat$subSystem_BC <- factor(combStat$subSystem_BC, levels = combStat$subSystem_BC)
	ggplot(combStat) +
		geom_point(aes(y = subSystem_BC, x = p.val.fdr),color='blue') +
		geom_vline(xintercept = 0.05, lty = 2) + theme_classic()+
		ggtitle("over-representation analysis, bd_lumped") +
		xlab("hypergeometric significance (fdr.adj.p.value)")	
	#write fdr.significant subset to csv	
	bd_lumped_subsystem_fdr = subset(combStat, p.val.fdr <= 0.05)
	write.table(bd_lumped_subsystem_fdr, "mta_tbl_prctile_top_significant_norm_t2/bd_lumped_subsystem_fdr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	bd_lumped_rxns_fdr = bd_lumped[bd_lumped$subSystem_BC %in% bd_lumped_subsystem_fdr$subSystem_BC,]
	write.table(bd_lumped_rxns_fdr, "mta_tbl_prctile_top_significant_norm_t2/bd_lumped_rxns_fdr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	#write all subsystems to csv
	BD_subsystem_all = combStat
	write.table(BD_subsystem_all, "mta_tbl_prctile_top_norm_t2/BD_subsystem_all.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	
# Fig.4: Number of rxns per subSystem (bd_responder)
	iPS_BD_R_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_R_TP_norm_t2.csv", header = T, sep = "\t")
	iPS_BD_R_TP_subSystem = iPS_BD_R_TP %>% count(SubSystem)
	bd_responder_subSystem = bd_responder %>% count(subSystem_BRC)
	comb = merge(bd_responder_subSystem, iPS_BD_R_TP_subSystem, by.x = "subSystem_BRC", by.y = "SubSystem")
	comb = comb[order(comb$n.x),]
	#4a
	dotchart(comb$n.x, labels = comb$subSystem_BRC,
			cex = 0.5, bg = "green", xlab = "# of rxns", main = "# of mta hits (bd_responder)")
	#4b
	dotchart(comb$n.y, labels = comb$subSystem_BRC,
			cex = 0.5, bg = "grey", xlab = "# of rxns", main = "# of mta hits (bd_responder) vs all rxns")	
	points(comb$n.x, 1:nrow(comb), col = "blue", cex = 0.7)	
	#4c
	p.val = phyper(q = comb$n.x-1, m = comb$n.y, n = nrow(iPS_BD_R_TP)-comb$n.y, k = sum(comb$n.x), lower.tail = FALSE, log.p = FALSE)
	p.val.fdr = p.adjust(p.val, method = "fdr")
	combStat = cbind(comb, p.val, p.val.fdr)
	combStat = combStat[order(-combStat$p.val.fdr),]
	combStat$subSystem_BRC <- factor(combStat$subSystem_BRC, levels = combStat$subSystem_BRC)
	ggplot(combStat) +
		geom_point(aes(y = subSystem_BRC, x = p.val.fdr),color='green') +
		geom_vline(xintercept = 0.05, lty = 2) + theme_classic()+
		ggtitle("over-representation analysis, bd_responder") +
		xlab("hypergeometric significance (fdr.adj.p.value)")
	#write fdr.significant subset to csv	
	bd_r_subsystem_fdr = subset(combStat, p.val.fdr <= 0.05)
	write.table(bd_r_subsystem_fdr, "mta_tbl_prctile_top_significant_norm_t2/bd_r_subsystem_fdr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	bd_r_rxns_fdr = bd_responder[bd_responder$subSystem_BRC %in% bd_r_subsystem_fdr$subSystem_BRC,]
	write.table(bd_r_rxns_fdr, "mta_tbl_prctile_top_significant_norm_t2/bd_r_rxns_fdr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	#write all subsystems to csv
	BD_R_subsystem_all = combStat
	write.table(BD_R_subsystem_all, "mta_tbl_prctile_top_norm_t2/BD_R_subsystem_all.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)

# Fig.5: Number of rxns per subSystem (bd_nonresponder)
	iPS_BD_NR_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_BD_NR_TP_norm_t2.csv", header = T, sep = "\t")
	iPS_BD_NR_TP_subSystem = iPS_BD_NR_TP %>% count(SubSystem)
	bd_nonresponder_subSystem = bd_nonresponder %>% count(subSystem_BNRC)
	comb = merge(bd_nonresponder_subSystem, iPS_BD_NR_TP_subSystem, by.x = "subSystem_BNRC", by.y = "SubSystem")
	comb = comb[order(comb$n.x),]
	#5a
	dotchart(comb$n.x, labels = comb$subSystem_BNRC,
			cex = 0.5, bg = "red", xlab = "# of rxns", main = "# of disrupted rxns (bd_nonresponder)")	
	#5b
	dotchart(comb$n.y, labels = comb$subSystem_BNRC,
			cex = 0.5, bg = "grey", xlab = "# of rxns", main = "# of disrupted rxns (bd_nonresponder) vs all rxns")	
	points(comb$n.x, 1:nrow(comb), col = "red", cex = 0.7)	
	#5c
	p.val = phyper(q = comb$n.x-1, m = comb$n.y, n = nrow(iPS_BD_NR_TP)-comb$n.y, k = sum(comb$n.x), lower.tail = FALSE, log.p = FALSE)
	p.val.fdr = p.adjust(p.val, method = "fdr")
	combStat = cbind(comb, p.val, p.val.fdr)
	combStat = combStat[order(-combStat$p.val.fdr),]
	combStat$subSystem_BNRC <- factor(combStat$subSystem_BNRC, levels = combStat$subSystem_BNRC)	
	ggplot(combStat) +
		geom_point(aes(y = subSystem_BNRC, x = p.val.fdr),color='red') +
		geom_vline(xintercept = 0.05, lty = 2) + theme_classic() +
		ggtitle("over-representation analysis, bd_nonresponder") +
		xlab("hypergeometric significance (fdr.adj.p.value)")	
	#write fdr.significant subset to csv	
	bd_nr_subsystem_fdr = subset(combStat, p.val.fdr <= 0.05)
	write.table(bd_nr_subsystem_fdr, "mta_tbl_prctile_top_significant_norm_t2/bd_nr_subsystem_fdr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	bd_nr_rxns_fdr = bd_nonresponder[bd_nonresponder$subSystem_BNRC %in% bd_nr_subsystem_fdr$subSystem_BNRC,]
	write.table(bd_nr_rxns_fdr, "mta_tbl_prctile_top_significant_norm_t2/bd_nr_rxns_fdr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	#write all subsystems to csv
	BD_NR_subsystem_all = combStat
	write.table(BD_NR_subsystem_all, "mta_tbl_prctile_top_norm_t2/BD_NR_subsystem_all.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)


# Fig.9: venn rxns (fdr.significant)
	x <- list(
		BD_Lumped = bd_lumped_rxns_fdr$del_rxnID_KO, 
		BD_Responder = bd_r_rxns_fdr$del_rxnID_KO, 
		BD_NonResponder = bd_nr_rxns_fdr$del_rxnID_KO)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlp in rxns (fdr.significant) disrupted between models")

# Fig.10: venn subSystems (fdr.significant)
	x <- list(
		BD_Lumped = bd_lumped_rxns_fdr$subSystem_BC, 
		BD_Responder = bd_r_rxns_fdr$subSystem_BRC, 
		BD_NonResponder = bd_nr_rxns_fdr$subSystem_BNRC)
	ggvenn(x, fill_color = c("blue", "green", "red")) + ggtitle("Overlp in subSystems (fdr.significant) disrupted between models")		

