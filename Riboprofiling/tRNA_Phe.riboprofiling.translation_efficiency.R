# 20210716 Phe brain translation efficiency analysis
####################
# TRANSLATION EFFICIENCY
####################

setwd("~/R/phe_translation_efficiency/")
library(stringr)
library(tidyr)
library(dplyr)
library(biomaRt)
library(ggplot2)

#===================
#### PREPARE DATA FILES ####
#===================

# import proteomic files and format
phe11_proteomics <- read.table("brain_phe11_directDIA_edited.txt", fill=T, sep="\t", quote="", header=T)
phe11_proteomics$AVG.Log2.Ratio <- -1 * phe11_proteomics$AVG.Log2.Ratio # switch FC as WT needs to be denominator

# import RNA files and format
phe11_rna <- read.table("brain_RNA_11all.txt", fill=T, sep="\t", quote="", header=T)
phe11_rna <- phe11_rna[,c(1:2,5:6,9)]

# join files 
phe11_merge <- left_join(phe11_proteomics, phe11_rna, by=c("Genes"="GeneName"))

# format
phe11_merge <- na.omit(phe11_merge)

#===================
#### ANALYSIS ####
#===================

# remove non-sig and calculate FC
phe11_merge <- phe11_merge[phe11_merge$Qvalue <= 0.05 | phe11_merge$padj <= 0.05,]
phe11_merge$pro.rna.Log2FC <- phe11_merge$AVG.Log2.Ratio - phe11_merge$log2FoldChange

# determine if FC has changed
phe11_merge$FC.value <- ifelse(phe11_merge$AVG.Log2.Ratio * phe11_merge$log2FoldChange < 0, "TRUE", "FALSE")


#===================
#### P FOR FINAL ####
#===================

# plot results highlighting those with changing FC only
FC.highlight.up <- phe11_merge[phe11_merge$pro.rna.Log2FC >=0.5,]
FC.highlight.down <- phe11_merge[phe11_merge$pro.rna.Log2FC <= -0.5,]
high.phe <- phe11_merge[grepl("Cds2|Syngr3|Synpr", phe11_merge$Genes),]

p3 <- phe11_merge %>% ggplot(aes(y=AVG.Log2.Ratio, x=log2FoldChange)) +
  geom_point(color="dark grey") +
  geom_point(data=FC.highlight.up, aes(y=AVG.Log2.Ratio, x=log2FoldChange), color="red") +
  geom_point(data=FC.highlight.down, aes(y=AVG.Log2.Ratio, x=log2FoldChange), color="blue") +
  labs(y="AVG.Log2.Ratio protein", x="log2FoldChange rna") +
  geom_text(aes(label=ifelse(Genes == "Calb2" | Genes == "Ncald" | Genes == "Slc32a1", as.character(Genes), ""), hjust=0, vjust=1)) +
  geom_abline(intercept=0.5, linetype="dashed") +
  geom_abline(intercept=-0.5, linetype="dashed") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), axis.line= element_line((color="black")))
p3 

ggsave(filename="efficiency_FC0.5.pdf", plot=p3, device="pdf")