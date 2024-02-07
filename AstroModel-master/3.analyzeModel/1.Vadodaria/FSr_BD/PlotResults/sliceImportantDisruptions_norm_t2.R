library("pheatmap")
library("RColorBrewer")
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
library(splitstackshape)
library(gplots)
library(ggplot2)
library(UpSetR)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(plot.matrix)

pdf("bd_tbl_norm_t2/FSR_UpSet_BD_norm_t2.pdf")

##########################################
# rnxs to matrix: import output
##########################################

	iPS_Ctrl_TP_vs_iPS_BD_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_TP_norm_t2.csv",header = T, sep = "\t");
	iPS_Ctrl_TP_vs_iPS_BD_R_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_R_TP_norm_t2.csv",header = T, sep = "\t");
	iPS_Ctrl_TP_vs_iPS_BD_NR_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/FSR_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_iPS_BD_NR_TP_norm_t2.csv",header = T, sep = "\t");
	UnChanged_iPS_Ctrl_TP_vs_Primary_TP = read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_Ctrl/UnChanged_iAstro_iPS_Ctrl_TP_norm_t2_vs_iAstro_Primary_TP_norm_t2.csv",header = T, sep = "\t");

	Ctrl_vs_BD = iPS_Ctrl_TP_vs_iPS_BD_TP$m
	Ctrl_vs_BD_R = iPS_Ctrl_TP_vs_iPS_BD_R_TP$m
	Ctrl_vs_BD_NR = iPS_Ctrl_TP_vs_iPS_BD_NR_TP$m
	Ctrl_vs_Primary_Unchanged = UnChanged_iPS_Ctrl_TP_vs_Primary_TP$n

###################################
# rnxs to matrix - wrangle
###################################

	mat = lst(Ctrl_vs_BD, Ctrl_vs_BD_R, Ctrl_vs_BD_NR, Ctrl_vs_Primary_Unchanged) %>% 
				
	  enframe %>% 
	  unnest %>% 
	  count(name, value) %>% 
	  spread(value, n, fill = 0)
	  
	mat = t(mat)

	write.table(mat, "bd_tbl_norm_t2/FSR_Table_BD_norm_t2.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

	mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/FSr_BD/PlotResults/bd_tbl_norm_t2/FSR_Table_BD_norm_t2.csv",
			header = T, sep = "\t")

##################################
# rnxs to matrix - UpSetR (TP)
##################################

	p = upset(mat,
		order.by = c("degree","freq"), #query.legend = "bottom",
		matrix.color="black", sets = c("Ctrl_vs_Primary_Unchanged", "Ctrl_vs_BD", "Ctrl_vs_BD_R", "Ctrl_vs_BD_NR"), 
		keep.order = TRUE, 
		queries = list(list(query = intersects, params = list("Ctrl_vs_BD", "Ctrl_vs_BD_R", "Ctrl_vs_BD_NR", "Ctrl_vs_Primary_Unchanged"), 
		color = "blue", active = T), 
		list(query = intersects, params = list("Ctrl_vs_BD", "Ctrl_vs_BD_R", "Ctrl_vs_Primary_Unchanged"), color = "blue", active = T),
		list(query = intersects, params = list("Ctrl_vs_BD", "Ctrl_vs_BD_NR", "Ctrl_vs_Primary_Unchanged"), color = "blue", active = T),
		list(query = intersects, params = list("Ctrl_vs_BD", "Ctrl_vs_Primary_Unchanged"), color = "blue", active = T),
		list(query = intersects, params = list("Ctrl_vs_BD_R", "Ctrl_vs_Primary_Unchanged"), color = "green", active = T),
		list(query = intersects, params = list("Ctrl_vs_BD_NR", "Ctrl_vs_Primary_Unchanged"), color = "red", active = T)),
		
		sets.bar.color=c("dark grey", "dark blue", "dark blue", "dark blue"), nintersects = 60)

	p

	svg(filename="bd_tbl_norm_t2/FSR_UpSet_BD_TP_norm_t2.svg")
	p 
	dev.off()

################
# rnxs slice
################

	#bd_lumped_a
	x = list(Ctrl_vs_Primary_Unchanged, Ctrl_vs_BD, Ctrl_vs_BD_R, Ctrl_vs_BD_NR)
	bd_lumped_a = Reduce(intersect, x)		
	length(bd_lumped_a)
	#~ write.table(bd_lumped_a, "bd_tbl_norm_t2/bd_lumped_a.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_lumped_b
	x = list(Ctrl_vs_Primary_Unchanged, Ctrl_vs_BD, Ctrl_vs_BD_R)
	xInt = Reduce(intersect, x)
	y = Ctrl_vs_BD_NR
	bd_lumped_b = setdiff(xInt,y)
	length(bd_lumped_b)
	#~ write.table(bd_lumped_b, "bd_tbl_norm_t2/bd_lumped_b.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_lumped_c
	x = list(Ctrl_vs_Primary_Unchanged, Ctrl_vs_BD, Ctrl_vs_BD_NR)
	xInt = Reduce(intersect, x)
	y = Ctrl_vs_BD_R
	bd_lumped_c = setdiff(xInt,y)
	length(bd_lumped_c)
	#~ write.table(bd_lumped_c, "bd_tbl_norm_t2/bd_lumped_c.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_lumped_d
	x = list(Ctrl_vs_Primary_Unchanged, Ctrl_vs_BD)
	xInt = Reduce(intersect,x)
	y = list(Ctrl_vs_BD_R, Ctrl_vs_BD_NR)
	yInt = Reduce(union,y)
	bd_lumped_d = setdiff(xInt,yInt)
	length(bd_lumped_d)
	#~ write.table(bd_lumped_d, "bd_tbl_norm_t2/bd_lumped_d.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_r
	x = list(Ctrl_vs_Primary_Unchanged, Ctrl_vs_BD_R)
	xInt = Reduce(intersect,x)
	y = list(Ctrl_vs_BD, Ctrl_vs_BD_NR)
	yInt = Reduce(union,y)
	bd_r = setdiff(xInt,yInt)
	length(bd_r)
	#~ write.table(bd_r, "bd_tbl_norm_t2/bd_r.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

	#bd_nr
	x = list(Ctrl_vs_Primary_Unchanged, Ctrl_vs_BD_NR)
	xInt = Reduce(intersect,x)
	y = list(Ctrl_vs_BD, Ctrl_vs_BD_R)
	yInt = Reduce(union,y)
	bd_nr = setdiff(xInt,yInt)
	length(bd_nr)
	#~ write.table(bd_nr, "bd_tbl_norm_t2/bd_nr.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)

#######################
# rnxs slice - subset
#######################

	#bd_lumped
	x = list(bd_lumped_a, bd_lumped_b, bd_lumped_c, bd_lumped_d)
	bd_lumped = Reduce(union,x)
	write.table(bd_lumped, "bd_tbl_norm_t2/bd_lumped.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)
	
	#bd_responder
	x = list(bd_r, bd_lumped_b)
	bd_responder = Reduce(union,x)
	write.table(bd_responder, "bd_tbl_norm_t2/bd_responder.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)
	
	#bd_nonresponder
	x = list(bd_nr, bd_lumped_c)
	bd_nonresponder = Reduce(union,x)
	write.table(bd_nonresponder, "bd_tbl_norm_t2/bd_nonresponder.csv", sep = "\t", quote = FALSE, row.names = F, col.names=F)
