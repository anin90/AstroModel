library("pheatmap")
library("RColorBrewer")
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
#install.packages("splitstackshape")
library(splitstackshape)
library(gplots)
library(ggplot2)
library(UpSetR)
library(grid)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
#install.packages("plot.matrix")
library(plot.matrix)

pdf("plot_fsr_intersect.pdf")

###############################################
#rnxs to matrix - import model constraints - v8
###############################################

Literature_Input = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/Literature_Input.csv");
Akkouh_v8_Input = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/Akkouh_v8_Input.csv");
GSE66276_Li_v8_Input = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/GSE66276_Li_v8_Input.csv");
GSE132397_Li_v8_Input = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/GSE132397_Li_v8_Input.csv");

Literature_Model_Input = Literature_Input
Rat_RNA_Seq_v8_Input = Akkouh_v8_Input
Rat_Microarray_v8_Input = GSE66276_Li_v8_Input
Mouse_Microarray_v8_Input = GSE132397_Li_v8_Input

##############################################
# rnxs to matrix: iMAT_model_TP_EXP_ASM_BBB_v8
##############################################

Literature_Out_v8 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/Literature_Out_v8.csv");
Akkouh_v8 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/Akkouh_v8.csv");
GSE66276_Li_v8 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/GSE66276_Li_v8.csv");
GSE132397_Li_v8 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/GSE132397_Li_v8.csv");

Literature_Model_v8_Output = Literature_Out_v8
Rat_RNA_Seq_v8_Output = Akkouh_v8
Rat_Microarray_v8_Output = GSE66276_Li_v8
Mouse_Microarray_v8_Output = GSE132397_Li_v8

###############################
#rnxs to matrix - wrangle - v8
###############################

mat_v8_rxns = lst(Literature_Model_v8_Output, Rat_RNA_Seq_v8_Output, Rat_Microarray_v8_Output, Mouse_Microarray_v8_Output,
					Literature_Model_Input, Rat_RNA_Seq_v8_Input, Rat_Microarray_v8_Input, Mouse_Microarray_v8_Input) %>% 
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)
  
mat_v8 = t(mat_v8_rxns)

write.table(mat_v8, "mat_v8.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

##############################
#rnxs to matrix - UpSetR - v8
##############################

mat_v8 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/mat_v8.csv",
		header = T, sep = "\t")

p = upset(mat_v8,
	order.by = c("degree","freq"), #query.legend = "bottom",
	matrix.color="black", sets = c("Literature_Model_v8_Output", "Rat_RNA_Seq_v8_Output", "Rat_Microarray_v8_Output", "Mouse_Microarray_v8_Output",
	"Literature_Model_Input", "Rat_RNA_Seq_v8_Input", "Rat_Microarray_v8_Input", "Mouse_Microarray_v8_Input"), keep.order = TRUE,
	sets.bar.color=c("dark grey","dark grey","dark grey","dark grey","light blue","light blue","light blue","light blue"),
	queries = list(list(query = intersects, 
	params = list("Rat_RNA_Seq_v8_Output", "Rat_Microarray_v8_Output", "Mouse_Microarray_v8_Output"), color = "maroon", active = T,
	query.name = "Output Intersections"),
	list(query = intersects, params = list("Mouse_Microarray_v8_Output", "Rat_Microarray_v8_Output"), color = "orange", active = T),
	list(query = intersects, params = list("Mouse_Microarray_v8_Output", "Rat_RNA_Seq_v8_Output"), color = "orange", active = T),
	list(query = intersects, params = list("Rat_Microarray_v8_Output", "Rat_RNA_Seq_v8_Output"), color = "orange", active = T),
	list(query = intersects, params = list("Rat_Microarray_v8_Output", "Mouse_Microarray_v8_Output", "Rat_RNA_Seq_v8_Input"), color = "orange", active = T),
	list(query = intersects, params = list("Rat_Microarray_v8_Input", "Mouse_Microarray_v8_Output", "Rat_RNA_Seq_v8_Output"), color = "orange", active = T)))
	
p 

###############################
#PLOT OUT ONLY
###############################

Literature = Literature_Model_v8_Output
Akkouh_et_al_Rat = Rat_RNA_Seq_v8_Output
GSE66276_Rat = Rat_Microarray_v8_Output
GSE132397_Mouse = Mouse_Microarray_v8_Output

###############################
#rnxs to matrix - wrangle - v8
###############################

mat_v8_rxns = lst(Literature, Akkouh_et_al_Rat, Rat_Microarray_v8_Output, GSE66276_Rat,
					GSE132397_Mouse) %>% 
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)
  
mat_v8 = t(mat_v8_rxns)

write.table(mat_v8, "mat_v8_out_only.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

##############################
#rnxs to matrix - UpSetR - v8
##############################

mat_v8 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/mat_v8_out_only.csv",
		header = T, sep = "\t")

p = upset(mat_v8,
	order.by = c("degree","freq"), #query.legend = "bottom",
	matrix.color="black", sets = c("Literature", "Akkouh_et_al_Rat", "GSE66276_Rat", "GSE132397_Mouse"), keep.order = TRUE,
	sets.bar.color=c("#F8CECC","#D5E8D4","#FFF2CC","#FFF2CC"),
	queries = list(list(query = intersects, 
	params = list("Akkouh_et_al_Rat", "GSE66276_Rat", "GSE132397_Mouse"), color = "maroon", active = T,
	query.name = "Output Intersections"),
	list(query = intersects, params = list("GSE132397_Mouse", "GSE66276_Rat"), color = "orange", active = T),
	list(query = intersects, params = list("GSE132397_Mouse", "Akkouh_et_al_Rat"), color = "orange", active = T),
	list(query = intersects, params = list("GSE66276_Rat", "Akkouh_et_al_Rat"), color = "orange", active = T)))
	
p 


############################
#intersection_heatmap - v8
############################

mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/mat_v8_intersection.csv",
		header = T, sep = "\t")
	
row.names(mat) <- mat$Rxn_Formula
mat <- mat[,3:5]
dat_matrix <- data.matrix(mat)
#col.pal <- brewer.pal(9,"Spectral")			#Spectral looks better
hm.parameters <- list(dat_matrix,cellwidth = 15, scale = "none",
cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 6, fontsize_col = 6,
color = c("light green", "white", "dark orange", "dark red"),breaks = c(0, 0.8, 1.5, 3, max(dat_matrix)))
do.call("pheatmap", hm.parameters)


mat <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/mat_v8_intersection.csv",
		header = T, sep = "\t")
	
row.names(mat) <- mat$Rxn_ID
mat <- mat[,3:5]
dat_matrix <- data.matrix(mat)
#col.pal <- brewer.pal(9,"Spectral")			#Spectral looks better
hm.parameters <- list(dat_matrix,cellwidth = 15, scale = "none",
cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 6, fontsize_col = 6,
color = c("light green", "white", "dark orange", "dark red"),breaks = c(0, 0.8, 1.5, 3, max(dat_matrix)))
do.call("pheatmap", hm.parameters)


###############################################
#rnxs to matrix - import model constraints - v9
###############################################

Literature_Input = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/Literature_Input.csv");
Akkouh_v9_Input = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/Akkouh_v9_Input.csv");
GSE66276_Li_v9_Input = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/GSE66276_Li_v9_Input.csv");
GSE132397_Li_v9_Input = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/GSE132397_Li_v9_Input.csv");

Literature_Model_Input = Literature_Input
Rat_RNA_Seq_v9_Input = Akkouh_v9_Input
Rat_Microarray_v9_Input = GSE66276_Li_v9_Input
Mouse_Microarray_v9_Input = GSE132397_Li_v9_Input

##############################################
# rnxs to matrix: iMAT_model_TP_EXP_ASM_BBB_v9
##############################################

Literature_Out_v9 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/Literature_Out_v9.csv");
Akkouh_v9 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/Akkouh_v9.csv");
GSE66276_Li_v9 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/GSE66276_Li_v9.csv");
GSE132397_Li_v9 = readLines("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/FSr_rxns/GSE132397_Li_v9.csv");

Literature_Model_v9_Output = Literature_Out_v9
Rat_RNA_Seq_v9_Output = Akkouh_v9
Rat_Microarray_v9_Output = GSE66276_Li_v9
Mouse_Microarray_v9_Output = GSE132397_Li_v9

###############################
#rnxs to matrix - wrangle - v9
###############################

mat_v9_rxns = lst(Literature_Model_v9_Output, Rat_RNA_Seq_v9_Output, Rat_Microarray_v9_Output, Mouse_Microarray_v9_Output,
					Literature_Model_Input, Rat_RNA_Seq_v9_Input, Rat_Microarray_v9_Input, Mouse_Microarray_v9_Input) %>% 
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)

mat_v9 = t(mat_v9_rxns)

write.table(mat_v9, "mat_v9.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=F)

############################
#rnxs to matrix - UpSetR - v9
############################

mat_v9 <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/Run/plot/mat_v9.csv",
		header = T, sep = "\t")

p = upset(mat_v9,
	order.by = c("degree","freq"),
	matrix.color="black", sets = c("Literature_Model_v9_Output", "Rat_RNA_Seq_v9_Output", "Rat_Microarray_v9_Output", "Mouse_Microarray_v9_Output",
	"Literature_Model_Input", "Rat_RNA_Seq_v9_Input", "Rat_Microarray_v9_Input", "Mouse_Microarray_v9_Input"), keep.order = TRUE,
	sets.bar.color=c("dark grey","dark grey","dark grey","dark grey","light blue","light blue","light blue","light blue"))
#	queries = list(list(query = intersects, 
#	params = list("Literature_curated_reaction", "Rat_Microarray_v9"), color = "orange", active = T),
#	list(query = intersects, params = list("Rat_Microarray_v9", "Rat_RNA_Seq_v9"), color = "orange", active = T)))

p 	
	

	
