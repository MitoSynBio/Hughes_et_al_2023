library(vsn)
library(apeglm)
library(tximport)
library(DESeq2)
library(EnhancedVolcano)

## Data
####################################################################################################################################
### Phe11 Brain ###
samplesPhe11Brain <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/brain.Phe11_vs_WT.n_4.samples", header=TRUE, sep="\t")
countsPhe11Brain <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/brain.Phe11_vs_WT.n_4.nfr.counts", header=TRUE, sep="\t")

ddsPhe11Brain <- DESeqDataSetFromMatrix(countsPhe11Brain,
                                          colData = samplesPhe11Brain,
                                          design = ~ Genotype)
ddsPhe11Brain$Genotype <- relevel(ddsPhe11Brain$Genotype, ref = "WT")
ddsPhe11Brain <- DESeq(ddsPhe11Brain)
resultsNames(ddsPhe11Brain)

### The effect of the Phe11 genotype - no threshold - s < 0.01 (svalues)
resPhe11_Brain <- results(ddsPhe11Brain, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe11_Brain <- lfcShrink(ddsPhe11Brain, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe11_Brain, svalue=TRUE)
write.table(resPhe11_Brain,"tables/Brain_Phe11_vs_WT.results.s_0.01.no_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11_Brain)

### The effect of the Phe11 genotype - log2(1.5) threshold - s < 0.01 (svalues)
resPhe11_Brain.thresh <- results(ddsPhe11Brain, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe11_Brain.thresh <- lfcShrink(ddsPhe11Brain, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe11_Brain.thresh, svalue=TRUE)
write.table(resPhe11_Brain.thresh,"tables/Brain_Phe11_vs_WT.results.s_0.01.log1.5_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11_Brain.thresh)

## Dispersion estimates
pdf("plots/Brain_tRNAPhe11_DispEstimates.indiv.nfr.pdf", width=10, height=10)
par(mfrow=c(1,1))
plotDispEsts(ddsPhe11Brain)
dev.off()

#ntdPhe11Brain <- normTransform(ddsPhe11Brain)
vsdPhe11Brain <- rlog(ddsPhe11Brain, blind=FALSE)
#rldPhe11Brain <- vst(ddsPhe11Brain, blind=FALSE)

#meanSdPlot(assay(ntdPhe11Brain))
#meanSdPlot(assay(vsdPhe11Brain))
#meanSdPlot(assay(rldPhe11Brain))

pdf("plots/Brain_tRNAPhe11_PCAplot.indiv.nfr.pdf", width=10, height=10)
plotPCA(vsdPhe11Brain, intgroup="Genotype")
dev.off()

# NormCounts
normCountsPhe11Brain <- counts(ddsPhe11Brain, normalized=TRUE)
write.table(normCountsPhe11Brain,"tables/Brain_tRNAPhe11_normCounts.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)



####################################################################################################################################
### Phe21 Brain ###
samplesPhe21Brain <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/brain.Phe21_vs_WT.n_3.samples", header=TRUE, sep="\t")
countsPhe21Brain <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/brain.Phe21_vs_WT.n_3.nfr.counts", header=TRUE, sep="\t")

ddsPhe21Brain <- DESeqDataSetFromMatrix(countsPhe21Brain,
                                          colData = samplesPhe21Brain,
                                          design = ~ Genotype)
ddsPhe21Brain$Genotype <- relevel(ddsPhe21Brain$Genotype, ref = "WT")
ddsPhe21Brain <- DESeq(ddsPhe21Brain)
resultsNames(ddsPhe21Brain)

### The effect of the Phe21 genotype - no threshold - s < 0.01 (svalues)
resPhe21_Brain <- results(ddsPhe21Brain, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe21_Brain <- lfcShrink(ddsPhe21Brain, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe21_Brain, svalue=TRUE)
write.table(resPhe21_Brain,"tables/Brain_Phe21_vs_WT.results.s_0.01.no_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21_Brain)

### The effect of the Phe21 genotype - log2(1.5) threshold - s < 0.01 (svalues)
resPhe21_Brain.thresh <- results(ddsPhe21Brain, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe21_Brain.thresh <- lfcShrink(ddsPhe21Brain, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe21_Brain.thresh, svalue=TRUE)
write.table(resPhe21_Brain.thresh,"tables/Brain_Phe21_vs_WT.results.s_0.01.log1.5_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21_Brain.thresh)

## Dispersion estimates
pdf("plots/Brain_tRNAPhe21_DispEstimates.indiv.nfr.pdf", width=10, height=10)
par(mfrow=c(1,1))
plotDispEsts(ddsPhe21Brain)
dev.off()

#ntdPhe21Brain <- normTransform(ddsPhe21Brain)
vsdPhe21Brain <- rlog(ddsPhe21Brain, blind=FALSE)
#rldPhe21Brain <- vst(ddsPhe21Brain, blind=FALSE)

#meanSdPlot(assay(ntdPhe21Brain))
#meanSdPlot(assay(vsdPhe21Brain))
#meanSdPlot(assay(rldPhe21Brain))

pdf("plots/Brain_tRNAPhe21_PCAplot.indiv.nfr.pdf", width=10, height=10)
plotPCA(vsdPhe21Brain, intgroup="Genotype")
dev.off()

# NormCounts
normCountsPhe21Brain <- counts(ddsPhe21Brain, normalized=TRUE)
write.table(normCountsPhe21Brain,"tables/Brain_tRNAPhe21_normCounts.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)

####################################################################################################################################

####################################################################################################################################
### Phe31 Brain ###
samplesPhe31Brain <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/brain.Phe31_vs_WT.n_3.samples", header=TRUE, sep="\t")
countsPhe31Brain <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/brain.Phe31_vs_WT.n_3.nfr.counts", header=TRUE, sep="\t")

ddsPhe31Brain <- DESeqDataSetFromMatrix(countsPhe31Brain,
                                          colData = samplesPhe31Brain,
                                          design = ~ Genotype)
ddsPhe31Brain$Genotype <- relevel(ddsPhe31Brain$Genotype, ref = "WT")
ddsPhe31Brain <- DESeq(ddsPhe31Brain)
resultsNames(ddsPhe31Brain)

### The effect of the Phe31 genotype - no threshold - s < 0.01 (svalues)
resPhe31_Brain <- results(ddsPhe31Brain, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe31_Brain <- lfcShrink(ddsPhe31Brain, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe31_Brain, svalue=TRUE)
write.table(resPhe31_Brain,"tables/Brain_Phe31_vs_WT.results.s_0.01.no_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe31_Brain)

### The effect of the Phe31 genotype - log1(1.5) threshold - s < 0.01 (svalues)
resPhe31_Brain.thresh <- results(ddsPhe31Brain, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe31_Brain.thresh <- lfcShrink(ddsPhe31Brain, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe31_Brain.thresh, svalue=TRUE)
write.table(resPhe31_Brain.thresh,"tables/Brain_Phe31_vs_WT.results.s_0.01.log1.5_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe31_Brain.thresh)

## Dispersion estimates
pdf("plots/Brain_tRNAPhe31_DispEstimates.indiv.nfr.pdf", width=10, height=10)
par(mfrow=c(1,1))
plotDispEsts(ddsPhe31Brain)
dev.off()

#ntdPhe31Brain <- normTransform(ddsPhe31Brain)
vsdPhe31Brain <- rlog(ddsPhe31Brain, blind=FALSE)
#rldPhe31Brain <- vst(ddsPhe31Brain, blind=FALSE)

#meanSdPlot(assay(ntdPhe31Brain))
#meanSdPlot(assay(vsdPhe31Brain))
#meanSdPlot(assay(rldPhe31Brain))

pdf("plots/Brain_tRNAPhe31_PCAplot.indiv.nfr.pdf", width=10, height=10)
plotPCA(vsdPhe31Brain, intgroup="Genotype")
dev.off()

# NormCounts
normCountsPhe31Brain <- counts(ddsPhe31Brain, normalized=TRUE)
write.table(normCountsPhe31Brain,"tables/Brain_tRNAPhe31_normCounts.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)

####################################################################################################################################

pdf("plots/Brain_tRNAPhe_MAplot.Phe_results.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11_Brain, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Brain, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe31_Brain, main="Phe31 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()


pdf("plots/Brain_tRNAPhe_MAplot.Phe_results.thresh.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11_Brain.thresh, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Brain.thresh, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe31_Brain.thresh, main="Phe31 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
save.image()


EnhancedVolcano(resPhe11_Brain,
                lab=rownames(resPhe11_Brain),
                x='log2FoldChange',
                y='svalue',
                pCutoff=0.01,
                FCcutoff=0,
                xlim=c(-2.5,2.5),
                drawConnectors=TRUE)

EnhancedVolcano(resPhe11_Brain.thresh,
                lab=rownames(resPhe11_Brain.thresh),
                x='log2FoldChange',
                y='svalue',
                pCutoff=0.01,
                FCcutoff=log2(1.5),
                xlim=c(-2.5,2.5),
                drawConnectors=TRUE)
