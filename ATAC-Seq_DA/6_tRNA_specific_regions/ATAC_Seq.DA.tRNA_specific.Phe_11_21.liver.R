library(vsn)
library(apeglm)
library(tximport)
library(DESeq2)

## Data
####################################################################################################################################
### Phe11 Liver ###
samplesPhe11Liver <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/liver.Phe11_vs_WT.n_4.samples", header=TRUE, sep="\t")
countsPhe11Liver <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/liver.Phe11_vs_WT.n_4.nfr.counts", header=TRUE, sep="\t")

ddsPhe11Liver <- DESeqDataSetFromMatrix(countsPhe11Liver,
                                        colData = samplesPhe11Liver,
                                        design = ~ Genotype)
ddsPhe11Liver$Genotype <- relevel(ddsPhe11Liver$Genotype, ref = "WT")
ddsPhe11Liver <- DESeq(ddsPhe11Liver)
resultsNames(ddsPhe11Liver)

### The effect of the Phe11 genotype - no threshold - s < 0.01 (svalues)
resPhe11_Liver <- results(ddsPhe11Liver, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe11_Liver <- lfcShrink(ddsPhe11Liver, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe11_Liver, svalue=TRUE)
write.table(resPhe11_Liver,"tables/Liver_Phe11_vs_WT.results.s_0.01.no_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11_Liver)

### The effect of the Phe11 genotype - log2(1.5) threshold - s < 0.01 (svalues)
resPhe11_Liver.thresh <- results(ddsPhe11Liver, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe11_Liver.thresh <- lfcShrink(ddsPhe11Liver, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe11_Liver.thresh, svalue=TRUE)
write.table(resPhe11_Liver.thresh,"tables/Liver_Phe11_vs_WT.results.s_0.01.log1.5_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11_Liver.thresh)

## Dispersion estimates
pdf("plots/Liver_tRNAPhe11_DispEstimates.indiv.nfr.pdf", width=10, height=10)
par(mfrow=c(1,1))
plotDispEsts(ddsPhe11Liver)
dev.off()

#ntdPhe11Liver <- normTransform(ddsPhe11Liver)
vsdPhe11Liver <- rlog(ddsPhe11Liver, blind=FALSE)
#rldPhe11Liver <- vst(ddsPhe11Liver, blind=FALSE)

#meanSdPlot(assay(ntdPhe11Liver))
#meanSdPlot(assay(vsdPhe11Liver))
#meanSdPlot(assay(rldPhe11Liver))

pdf("plots/Liver_tRNAPhe11_PCAplot.indiv.nfr.pdf", width=10, height=10)
plotPCA(vsdPhe11Liver, intgroup="Genotype")
dev.off()

# NormCounts
normCountsPhe11Liver <- counts(ddsPhe11Liver, normalized=TRUE)
write.table(normCountsPhe11Liver,"tables/Liver_tRNAPhe11_normCounts.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)



####################################################################################################################################
### Phe21 Liver ###
samplesPhe21Liver <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/liver.Phe21_vs_WT.n_3.samples", header=TRUE, sep="\t")
countsPhe21Liver <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/liver.Phe21_vs_WT.n_3.nfr.counts", header=TRUE, sep="\t")

ddsPhe21Liver <- DESeqDataSetFromMatrix(countsPhe21Liver,
                                        colData = samplesPhe21Liver,
                                        design = ~ Genotype)
ddsPhe21Liver$Genotype <- relevel(ddsPhe21Liver$Genotype, ref = "WT")
ddsPhe21Liver <- DESeq(ddsPhe21Liver)
resultsNames(ddsPhe21Liver)

### The effect of the Phe21 genotype - no threshold - s < 0.01 (svalues)
resPhe21_Liver <- results(ddsPhe21Liver, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe21_Liver <- lfcShrink(ddsPhe21Liver, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe21_Liver, svalue=TRUE)
write.table(resPhe21_Liver,"tables/Liver_Phe21_vs_WT.results.s_0.01.no_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21_Liver)

### The effect of the Phe21 genotype - log2(1.5) threshold - s < 0.01 (svalues)
resPhe21_Liver.thresh <- results(ddsPhe21Liver, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe21_Liver.thresh <- lfcShrink(ddsPhe21Liver, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe21_Liver.thresh, svalue=TRUE)
write.table(resPhe21_Liver.thresh,"tables/Liver_Phe21_vs_WT.results.s_0.01.log1.5_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21_Liver.thresh)

## Dispersion estimates
pdf("plots/Liver_tRNAPhe21_DispEstimates.indiv.nfr.pdf", width=10, height=10)
par(mfrow=c(1,1))
plotDispEsts(ddsPhe21Liver)
dev.off()

#ntdPhe21Liver <- normTransform(ddsPhe21Liver)
vsdPhe21Liver <- rlog(ddsPhe21Liver, blind=FALSE)
#rldPhe21Liver <- vst(ddsPhe21Liver, blind=FALSE)

#meanSdPlot(assay(ntdPhe21Liver))
#meanSdPlot(assay(vsdPhe21Liver))
#meanSdPlot(assay(rldPhe21Liver))

pdf("plots/Liver_tRNAPhe21_PCAplot.indiv.nfr.pdf", width=10, height=10)
plotPCA(vsdPhe21Liver, intgroup="Genotype")
dev.off()

# NormCounts
normCountsPhe21Liver <- counts(ddsPhe21Liver, normalized=TRUE)
write.table(normCountsPhe21Liver,"tables/Liver_tRNAPhe21_normCounts.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)

####################################################################################################################################

####################################################################################################################################

####################################################################################################################################

pdf("plots/Liver_tRNAPhe_MAplot.Phe_results.pdf", width=15, height=10)
par(mfrow=c(1,2))
plotMA(resPhe11_Liver, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Liver, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()


pdf("plots/Liver_tRNAPhe_MAplot.Phe_results.thresh.pdf", width=15, height=10)
par(mfrow=c(1,2))
plotMA(resPhe11_Liver.thresh, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Liver.thresh, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
save.image()



