#~ install.packages("pacman")
#~ install.packages("mgsub")
#~ install.packages("gsubfn")
#~ install.packages("readxl")
#~ install.packages("openxlsx")
pacman::p_load(pheatmap, RColorBrewer, dplyr, tidyverse, tidyr, tibble, splitstackshape, gplots, ggplot2, ggfortify, reshape2, 
factoextra, plot.matrix, VennDiagram, ggvenn, plotrix, pheatmap, magrittr, venn, mgsub, gsubfn, readxl, openxlsx, UpSetR) 

pdf("mta_tbl_prctile_top_abs/analyzeMTAscores_abs.pdf")
  
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
#~ Load koskuvi_mta
#######################

# Load data
	Ctrl_to_Primary = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_abs/Koskuvi_Ctrl_to_Primary.csv", header = T, sep = "\t")
	Primary_to_Ctrl = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_abs/Koskuvi_Primary_to_Ctrl.csv", header = T, sep = "\t")
	Ctrl_to_HT = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_abs/Koskuvi_Ctrl_to_HT.csv", header = T, sep = "\t")
	HT_to_Ctrl = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_abs/Koskuvi_HT_to_Ctrl.csv", header = T, sep = "\t")
	Ctrl_to_ST = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_abs/Koskuvi_Ctrl_to_ST.csv", header = T, sep = "\t")
	ST_to_Ctrl = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_abs/Koskuvi_ST_to_Ctrl.csv", header = T, sep = "\t")

# Rename data
	CP = Ctrl_to_Primary
	PC = Primary_to_Ctrl
	CH = Ctrl_to_HT
	HC = HT_to_Ctrl
	CS = Ctrl_to_ST
	SC = ST_to_Ctrl

# Histogram mta_score
	par(mfrow=c(4,2))
	hist(CP$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Ctrl_to_Primary")
	hist(PC$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Primary_to_Ctrl")
	hist(CH$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Ctrl_to_HT")
	hist(HC$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="HT_to_Ctrl")
	hist(CS$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="Ctrl_to_ST")
	hist(SC$mta_score, xlab = "score", ylab = "freq", col=rgb(1,0,0,0.5), border = "blue", breaks=30, main="ST_to_Ctrl")

# FilterMTAscore, all
	FilterMTAscore <- function(x) {
		mta_score_prctile <- ecdf(x$mta_score)(x$mta_score)
		mta_score_prctile <- mta_score_prctile*100
		x <- cbind(x, mta_score_prctile)
		x <- arrange(x, -mta_score_prctile)
		return(x) 
		}		

	dat = list(CP, PC, CH, HC, CS, SC)
	dat_all = lapply(dat, FilterMTAscore)
	
	CP_prctile = dat_all[[1]]
	PC_prctile = dat_all[[2]]
	CH_prctile = dat_all[[3]]
	HC_prctile = dat_all[[4]]
	CS_prctile = dat_all[[5]]
	SC_prctile = dat_all[[6]]
	
	write.table(dat_all[[1]], "mta_tbl_prctile_all_abs/Ctrl_to_Primary.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[2]], "mta_tbl_prctile_all_abs/Primary_to_Ctrl.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[3]], "mta_tbl_prctile_all_abs/Ctrl_to_HT.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[4]], "mta_tbl_prctile_all_abs/HT_to_Ctrl.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[5]], "mta_tbl_prctile_all_abs/Ctrl_to_ST.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	write.table(dat_all[[6]], "mta_tbl_prctile_all_abs/ST_to_Ctrl.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)


# FilterMTAscore, bottom 80% (unchanged)
	FilterMTAscore_unchanged <- function(x) {
		mta_score_prctile <- ecdf(x$mta_score)(x$mta_score)
		mta_score_prctile <- mta_score_prctile*100
		x <- cbind(x, mta_score_prctile)
		x <- arrange(x, -mta_score_prctile)
		x <- filter(x, mta_score_prctile < 80)
		return(x) 
		}		
	dat = list(CP, PC, CH, HC)
	dat_unchanged = lapply(dat, FilterMTAscore_unchanged)

	CP_filt = dat_unchanged[[1]]
	PC_filt = dat_unchanged[[2]]
	CH_filt = dat_unchanged[[3]]
	HC_filt = dat_unchanged[[4]]

# FilterMTAscore, top 20% (significant)
	FilterMTAscore_significant <- function(x) {
		mta_score_prctile <- ecdf(x$mta_score)(x$mta_score)
		mta_score_prctile <- mta_score_prctile*100
		x <- cbind(x, mta_score_prctile)
		x <- arrange(x, -mta_score_prctile)
		x <- filter(x, mta_score_prctile >= 80)
		return(x) 
		}	
	dat = list(CS, SC)			
	dat_significant = lapply(dat, FilterMTAscore_significant)

	CS_filt = dat_significant[[1]]
	SC_filt = dat_significant[[2]]
	
##################
# Load MTA rxns
##################

	CP_rxns = CP_filt$del_rxnID_KO
	PC_rxns = PC_filt$del_rxnID_KO
	CH_rxns = CH_filt$del_rxnID_KO
	HC_rxns = HC_filt$del_rxnID_KO
	CS_rxns = CS_filt$del_rxnID_KO
	SC_rxns = SC_filt$del_rxnID_KO

################
# rnxs slice
################

	Ctrl_unchanged = Reduce(intersect, list(CP_rxns, PC_rxns, CH_rxns, HC_rxns))
	ST_top = union(CS_rxns, SC_rxns)

###################################
# rnxs to matrix - wrangle
###################################

	mat = lst(Ctrl_unchanged, ST_top) %>% 
				
	  enframe %>% 
	  unnest %>% 
	  count(name, value) %>% 
	  spread(value, n, fill = 0)

	mat = t(mat)

	write.table(mat, "mta_tbl_prctile_top_abs/MTA_Table_ST_abs.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

	mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/PlotResults/mta_tbl_prctile_top_abs/MTA_Table_ST_abs.csv",
			header = T, sep = "\t")

##################################
# rnxs to matrix - UpSetR (TP)
##################################

	p = upset(mat,
		order.by = c("degree","freq"), #query.legend = "bottom",
		matrix.color="black", sets = c("Ctrl_unchanged", "ST_top"), 
		keep.order = TRUE, 
		queries = list(list(query = intersects,  params = list("Ctrl_unchanged", "ST_top"), 
		color = "red", active = F)), sets.bar.color=c("dark grey", "dark blue"), nintersects = 60)

	p

################
# rnxs slice
################

	#st_rxns
	x = list(Ctrl_unchanged, ST_top)
	st_rxns = Reduce(intersect, x)		
	length(st_rxns)
	write.table(st_rxns, "mta_tbl_prctile_top_abs/st_rxns.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)
	
################################
# grep prctile for sliced rxns
################################

	#all_rxns
	df_list = list(CP_prctile, PC_prctile, CH_prctile, HC_prctile, CS_prctile, SC_prctile)
					
	df_list = df_list %>% reduce(full_join, by='del_rxnID_KO')
	df_list_prctile = df_list[,grepl("*mta_score_prctile",names(df_list))]
	df_list_subSystem = df_list[,grepl("*subSystem",names(df_list))]
	df_list_MetabolicUnits = df_list[,grepl("*MetabolicUnits",names(df_list))]
	df_list_Localization = df_list[,grepl("*Localization",names(df_list))]
	df_list_subset = cbind(df_list$del_rxnID_KO, df_list_prctile, 
							df_list_subSystem, df_list_MetabolicUnits, df_list_Localization)
	
	colnames(df_list_subset) = c("del_rxnID_KO", "prctile_CP", "prctile_PC", "prctile_CH", "prctile_HC", "prctile_CS", "prctile_SC",
					"subSystem_CP", "subSystem_PC", "subSystem_CH", "subSystem_HC", "subSystem_CS", "subSystem_SC",
					"MetabolicUnits_CP", "MetabolicUnits_PC", "MetabolicUnits_CH", "MetabolicUnits_HC", "MetabolicUnits_CS", "MetabolicUnits_SC",
					"Localization_CP", "Localization_PC", "Localization_CH", "Localization_HC", "Localization_CS", "Localization_SC")
					
	prctile_tbl_all = df_list_subset

	#st_rxns
	st_rxns = prctile_tbl_all[prctile_tbl_all$del_rxnID_KO %in% st_rxns ,]
	write.table(st_rxns, "mta_tbl_prctile_top_abs/st_rxns.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)

###########
# Plot data
###########
			
# Fig.1: Number of rxns per subSystem (st_rxns)
	par(mfrow=c(1,1))
	iPS_ST_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Annotations/AnnotateRxnSubsystems/Rxns_iPS_SCZ_ST_TP_abs.csv", header = T, sep = "\t")
	iPS_ST_TP_subSystem = iPS_ST_TP %>% count(SubSystem)
	ST_subSystem = st_rxns %>% count(subSystem_SC)
	comb = merge(ST_subSystem, iPS_ST_TP_subSystem, by.x = "subSystem_SC", by.y = "SubSystem")
	comb = comb[order(comb$n.x),]
	#1a
	dotchart(comb$n.x, labels = comb$subSystem_SC,
			cex = 0.5, bg = "blue", xlab = "# of rxns", main = "# of mta hits (st_rxns)")
	#1b
	dotchart(comb$n.y, labels = comb$subSystem_SC,
			cex = 0.5, bg = "grey", xlab = "# of rxns", main = "# of mta hits (st_rxns) vs all rxns")	
	points(comb$n.x, 1:nrow(comb), col = "blue", cex = 0.7)	
	#1c
	p.val = phyper(q = comb$n.x-1, m = comb$n.y, n = nrow(iPS_ST_TP)-comb$n.y, k = sum(comb$n.x), lower.tail = FALSE, log.p = FALSE)
	p.val.fdr = p.adjust(p.val, method = "fdr")
	combStat = cbind(comb, p.val, p.val.fdr)
	combStat = combStat[order(-combStat$p.val.fdr),]
	combStat$subSystem_SC <- factor(combStat$subSystem_SC, levels = combStat$subSystem_SC)
	ggplot(combStat) +
		geom_point(aes(y = subSystem_SC, x = p.val.fdr),color='blue') +
		geom_vline(xintercept = 0.05, lty = 2) + theme_classic()+
		ggtitle("over-representation analysis, bd_lumped") +
		xlab("hypergeometric significance (fdr.adj.p.value)")
	#write fdr.significant subset to csv	
	ST_subsystem_fdr = subset(combStat, p.val.fdr <= 0.05)
	write.table(ST_subsystem_fdr, "mta_tbl_prctile_top_significant_abs/ST_subsystem_fdr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	ST_rxns_fdr = st_rxns[st_rxns$subSystem_SC %in% ST_subsystem_fdr$subSystem_SC,]
	write.table(ST_rxns_fdr, "mta_tbl_prctile_top_significant_abs/ST_rxns_fdr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)
	#write all subsystems to csv
	ST_subsystem_all = combStat
	write.table(ST_subsystem_all, "mta_tbl_prctile_top_abs/ST_subsystem_all.csv", sep = "\t", quote = FALSE, row.names = F, col.names=T)

