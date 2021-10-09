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
library(scales)
library(plyr)
pdf("Test4HumanAstro.pdf")

###########
# Load data
###########

Test4HumanAstro_TP <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Test4HumanAstro_TP.csv",
		header = T, sep = "\t")
Test4HumanAstro_T <- read.csv("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_T_TP/Test4HumanAstro_T.csv",
		header = T, sep = "\t")

# ATP Synthesis
#TP
ss.data = Test4HumanAstro_TP[Test4HumanAstro_TP$TestSolutionGroup == 'ATP demand',]
mdat = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("Primary_TP", "iPS_Ctrl_TP", "iPS_BD_TP", "iPS_BD_R_TP", "iPS_BD_NR_TP"))
ggplot(mdat, aes(x = variable, y = value))+
  geom_col(aes(fill = TestedMetabolite), width = 0.7, colour = "black") + theme_classic()+ 
   labs(fill = "Nutrient source", x ="closedModel", y = "ATP demand (mmol/gDw/h), per 1 mmol of nutrient uptake") + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_brewer(palette = "Pastel1")

ggplot(mdat, aes(x = variable, y = value, fill = TestedMetabolite)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1") + theme_classic() + 
   labs(fill = "Nutrient source", x ="closedModel", y = "ATP demand (mmol/gDw/h), per 1 mmol of nutrient uptake") + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#T
ss.data = Test4HumanAstro_T[Test4HumanAstro_T$TestSolutionGroup == 'ATP demand',]
mdat = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("Primary_T", "iPS_Ctrl_T", "iPS_BD_T", "iPS_BD_R_T", "iPS_BD_NR_T"))
ggplot(mdat, aes(x = variable, y = value))+
  geom_col(aes(fill = TestedMetabolite), width = 0.7, colour = "black") + theme_classic()+ 
   labs(fill = "Nutrient source", x ="closedModel", y = "ATP demand (mmol/gDw/h), per 1 mmol of nutrient uptake") + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_brewer(palette = "Pastel1")

ggplot(mdat, aes(x = variable, y = value, fill = TestedMetabolite)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1") + theme_classic() + 
   labs(fill = "Nutrient source", x ="closedModel", y = "ATP demand (mmol/gDw/h), per 1 mmol of nutrient uptake") + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))   
   

# Pyruvate utilization
#TP
ss.data = Test4HumanAstro_TP[Test4HumanAstro_TP$TestSolutionGroup == 'Pyruvate utilization',]
mdat = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("Primary_TP", "iPS_Ctrl_TP", "iPS_BD_TP", "iPS_BD_R_TP", "iPS_BD_NR_TP"))
ggplot(mdat, aes(x = variable, y = value))+
  geom_col(aes(fill = TestedMetabolite), width = 0.7, colour = "black") + theme_classic()+ 
   labs(fill = "Nutrient synthesis", x ="closedModel", y = "Pyruvate utilization (mmol/gDw/h), per 1 mmol of nutrient uptake") + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_brewer(palette = "Pastel1")

ggplot(mdat, aes(x = variable, y = value, fill = TestedMetabolite)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1") + theme_classic() + 
   labs(fill = "Nutrient synthesis", x ="closedModel", y = "Pyruvate utilization (mmol/gDw/h), per 1 mmol of nutrient uptake") + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#T
ss.data = Test4HumanAstro_T[Test4HumanAstro_T$TestSolutionGroup == 'Pyruvate utilization',]
mdat = melt(ss.data, id.vars=c("TestedMetabolite", "TestSolutionGroup"),
            measure.vars=c("Primary_T", "iPS_Ctrl_T", "iPS_BD_T", "iPS_BD_R_T", "iPS_BD_NR_T"))
ggplot(mdat, aes(x = variable, y = value))+
  geom_col(aes(fill = TestedMetabolite), width = 0.7, colour = "black") + theme_classic()+ 
   labs(fill = "Nutrient synthesis", x ="closedModel", y = "Pyruvate utilization (mmol/gDw/h), per 1 mmol of nutrient uptake") + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_brewer(palette = "Pastel1")

ggplot(mdat, aes(x = variable, y = value, fill = TestedMetabolite)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1") + theme_classic() + 
   labs(fill = "Nutrient synthesis", x ="closedModel", y = "Pyruvate utilization (mmol/gDw/h), per 1 mmol of nutrient uptake") + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






