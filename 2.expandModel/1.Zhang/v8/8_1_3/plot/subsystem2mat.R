library("pheatmap")
library("RColorBrewer")
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
install.packages("splitstackshape")
library(splitstackshape)

# Read the data into R
#df = read.table("/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/2_expandModel/8_1_3/plot/Li_subsystem.txt", sep="\t", header=T);

MI1PP=c("Inositol phosphate metabolism")
INSTt4=c()
PGS=c()
PDE1=c("Nucleotide interconversion")
BPNT=c()

PGMT=c("Sphingolipid metabolism",
"Starch and sucrose metabolism",
"Fructose and mannose metabolism",
"Transport, extracellular",
"Glycolysis/gluconeogenesis")

ADNCYC=c("Nucleotide interconversion")
MI1PP_INSTt4=c("Inositol phosphate metabolism")

Literature=c("Sphingolipid metabolism",
"Starch and sucrose metabolism",
"Fructose and mannose metabolism",
"Inositol phosphate metabolism",
"Transport, extracellular",
"Glycolysis/gluconeogenesis")

T13_F3=c("Folate metabolism",
"Nucleotide interconversion",
"Transport, peroxisomal",
"R group synthesis",
"Exchange/demand reaction",
"Butanoate metabolism",
"Transport, extracellular",
"Valine, leucine, and isoleucine metabolism",
"Transport, mitochondrial",
"Glycine, serine, alanine, and threonine metabolism",
"Transport, lysosomal",
"Miscellaneous",
"Fatty acid oxidation",
"Triacylglycerol synthesis",
"Transport, endoplasmic reticular",
"Fatty acid synthesis")

T13_F2=c("Folate metabolism",
"Glycerophospholipid metabolism",
"Arachidonic acid metabolism",
"Nucleotide interconversion",
"Transport, peroxisomal",
"R group synthesis",
"Exchange/demand reaction",
"Glyoxylate and dicarboxylate metabolism",
"Butanoate metabolism",
"Transport, extracellular",
"Eicosanoid metabolism",
"Valine, leucine, and isoleucine metabolism",
"Transport, mitochondrial",
"Glycine, serine, alanine, and threonine metabolism",
"Transport, lysosomal",
"Miscellaneous",
"Urea cycle",
"Fatty acid synthesis",
"Transport, endoplasmic reticular",
"Fatty acid oxidation",
"Triacylglycerol synthesis",
"Leukotriene metabolism")

T14_F3=c("Glycerophospholipid metabolism",
"Transport, golgi apparatus",
"Arachidonic acid metabolism",
"Eicosanoid metabolism",
"Transport, mitochondrial",
"Cholesterol metabolism",
"Bile acid synthesis",
"Nucleotide interconversion",
"Exchange/demand reaction",
"Glyoxylate and dicarboxylate metabolism",
"Fatty acid oxidation",
"Leukotriene metabolism",
"Transport, extracellular",
"Methionine and cysteine metabolism")

T14_F2=c("Glycerophospholipid metabolism",
"Arachidonic acid metabolism",
"Citric acid cycle",
"Transport, peroxisomal",
"Nucleotide interconversion",
"Exchange/demand reaction",
"Glyoxylate and dicarboxylate metabolism",
"Oxidative phosphorylation",
"Inositol phosphate metabolism",
"Transport, extracellular",
"Purine synthesis",
"Transport, golgi apparatus",
"ROS detoxification",
"Nucleotide salvage pathway",
"Eicosanoid metabolism",
"Transport, mitochondrial",
"Valine, leucine, and isoleucine metabolism",
"Taurine and hypotaurine metabolism",
"Bile acid synthesis",
"Transport, nuclear",
"Miscellaneous",
"Pyruvate metabolism",
"Fatty acid synthesis",
"Transport, endoplasmic reticular",
"Fatty acid oxidation",
"Leukotriene metabolism",
"Methionine and cysteine metabolism")

GSE66276_Li_F3=c("Purine synthesis",
"Exchange/demand reaction",
"Fatty acid oxidation",
"Inositol phosphate metabolism",
"Nucleotide interconversion",
"Transport, extracellular",
"Methionine and cysteine metabolism")

GSE66276_Li_F2 = c("Transport, peroxisomal",
"Nucleotide interconversion",
"Exchange/demand reaction",
"Inositol phosphate metabolism",
"Transport, extracellular",
"Purine synthesis",
"Glutathione metabolism",
"Galactose metabolism",
"Transport, mitochondrial",
"Transport, lysosomal",
"Sphingolipid metabolism",
"Phosphatidylinositol phosphate metabolism",
"Transport, endoplasmic reticular",
"Fatty acid oxidation",
"Fatty acid synthesis",
"Methionine and cysteine metabolism")

GSE132397_Li_F3 = c("Urea cycle",
"Arginine and proline metabolism",
"Exchange/demand reaction")

GSE132397_Li_F2 = c("Urea cycle",
"Arginine and proline metabolism",
"Exchange/demand reaction",
"Inositol phosphate metabolism",
"Transport, extracellular")

mat = lst(Literature, T13_F3, T13_F2, GSE66276_Li_F3, GSE66276_Li_F2, GSE132397_Li_F3, GSE132397_Li_F2) %>% 
  enframe %>% 
  unnest %>% 
  count(name, value) %>% 
  spread(value, n, fill = 0)
  
mat = t(mat)

write.table(mat, "subsystem2mat.csv", sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
