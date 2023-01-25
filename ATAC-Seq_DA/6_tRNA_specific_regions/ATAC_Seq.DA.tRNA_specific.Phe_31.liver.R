library(vsn)
library(apeglm)
library(tximport)
library(DESeq2)

## Data
####################################################################################################################################
### Phe31 Liver ###
samplesPhe31Liver <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/liver.Phe31_vs_WT.n_3.samples", header=TRUE, sep="\t")
countsPhe31Liver <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/7_tRNA_specific/liver.Phe31_vs_WT.n_3.nfr.counts", header=TRUE, sep="\t")

ddsPhe31Liver <- DESeqDataSetFromMatrix(countsPhe31Liver,
                                        colData = samplesPhe31Liver,
                                        design = ~ Genotype)
ddsPhe31Liver$Genotype <- relevel(ddsPhe31Liver$Genotype, ref = "WT")
ddsPhe31Liver <- DESeq(ddsPhe31Liver)
resultsNames(ddsPhe31Liver)

### The effect of the Phe31 genotype - no threshold - s < 0.01 (svalues)
resPhe31_Liver <- results(ddsPhe31Liver, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe31_Liver <- lfcShrink(ddsPhe31Liver, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe31_Liver, svalue=TRUE)
write.table(resPhe31_Liver,"tables/Liver_Phe31_vs_WT.results.s_0.01.no_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe31_Liver)

### The effect of the Phe31 genotype - log2(1.5) threshold - s < 0.01 (svalues)
resPhe31_Liver.thresh <- results(ddsPhe31Liver, name="Genotype_KO_vs_WT", alpha=0.01)
resPhe31_Liver.thresh <- lfcShrink(ddsPhe31Liver, coef="Genotype_KO_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe31_Liver.thresh, svalue=TRUE)
write.table(resPhe31_Liver.thresh,"tables/Liver_Phe31_vs_WT.results.s_0.01.log1.5_thresh.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe31_Liver.thresh)

## Dispersion estimates
pdf("plots/Liver_tRNAPhe31_DispEstimates.indiv.nfr.pdf", width=10, height=10)
par(mfrow=c(1,1))
plotDispEsts(ddsPhe31Liver)
dev.off()

#ntdPhe31Liver <- normTransform(ddsPhe31Liver)
vsdPhe31Liver <- rlog(ddsPhe31Liver, blind=FALSE)
#rldPhe31Liver <- vst(ddsPhe31Liver, blind=FALSE)

#meanSdPlot(assay(ntdPhe31Liver))
#meanSdPlot(assay(vsdPhe31Liver))
#meanSdPlot(assay(rldPhe31Liver))

pdf("plots/Liver_tRNAPhe31_PCAplot.indiv.nfr.pdf", width=10, height=10)
plotPCA(vsdPhe31Liver, intgroup="Genotype")
dev.off()

# NormCounts
normCountsPhe31Liver <- counts(ddsPhe31Liver, normalized=TRUE)
write.table(normCountsPhe31Liver,"tables/Liver_tRNAPhe31_normCounts.indiv.nfr.txt", sep="\t", col.names=T, row.names=T, quote=F)

####################################################################################################################################

####################################################################################################################################

pdf("plots/Liver_tRNAPhe_MAplot.Phe_results.updated.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11_Liver, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Liver, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe31_Liver, main="Phe31 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()


pdf("plots/Liver_tRNAPhe_sMAplot.Phe_results.thresh.updated.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11_Liver.thresh, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Liver.thresh, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe31_Liver.thresh, main="Phe31 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
save.image()

