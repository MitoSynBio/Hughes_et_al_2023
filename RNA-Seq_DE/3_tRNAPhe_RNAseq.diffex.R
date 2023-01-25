library(vsn)
library(DESeq2)
library(apeglm)
library(tximport)

## Data
dir <- "/media/Data/s.siira/tRNA_Phe_RNAseq/pipeline/2_alignment"
samples <- read.table("/media/Data/s.siira/tRNA_Phe_RNAseq/pipeline/2_alignment/tRNAPhe_samples.txt", header=TRUE)
files <- file.path(dir, samples$Location, "quant.sf")
names(files) <- samples$Sample

tx2gene <- read.table("/media/Data/s.siira/tRNA_Phe_RNAseq/pipeline/3_diffEx/t2g.txt", header=T, sep="\t", quote='')


########################################################################################################################################################
### Liver ###
########################################################################################################################################################
samplesLiver <- subset(samples, Tissue=="Liver")
txiLiver <- tximport(files[1:12], type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")

ddsLiver <- DESeqDataSetFromTximport(txiLiver,
                                colData = samplesLiver,
                                design = ~ Genotype)
ddsLiver$Genotype <- relevel(ddsLiver$Genotype, ref = "WT")
ddsLiver <- DESeq(ddsLiver)
resultsNames(ddsLiver)

### The effect of the Phe11 genotype - no threshold - s < 0.01 (svalues)
resPhe11_Liver <- results(ddsLiver, name="Genotype_Phe11_vs_WT", alpha=0.01)
resPhe11_Liver <- lfcShrink(ddsLiver, coef="Genotype_Phe11_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe11_Liver, svalue=TRUE)
write.table(resPhe11_Liver,"tables/Liver_Phe11_vs_WT.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11_Liver)

resPhe21_Liver <- results(ddsLiver, name="Genotype_Phe21_vs_WT", alpha=0.01)
resPhe21_Liver <- lfcShrink(ddsLiver, coef="Genotype_Phe21_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe21_Liver, svalue=TRUE)
write.table(resPhe21_Liver,"tables/Liver_Phe21_vs_WT.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21_Liver)

resPhe31_Liver <- results(ddsLiver, name="Genotype_Phe31_vs_WT", alpha=0.01)
resPhe31_Liver <- lfcShrink(ddsLiver, coef="Genotype_Phe31_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe31_Liver, svalue=TRUE)
write.table(resPhe31_Liver,"tables/Liver_Phe31_vs_WT.results.s0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe31_Liver)


## Dispersion estimates
pdf("plots/Liver_tRNAPhe_DispEstimates.pdf", width=10, height=10)
plotDispEsts(ddsLiver)
dev.off()

ntdLiver <- normTransform(ddsLiver)
rldLiver <- rlog(ddsLiver, blind=FALSE)
vsdLiver <- vst(ddsLiver, blind=FALSE)

meanSdPlot(assay(ntdLiver))
meanSdPlot(assay(vsdLiver))
meanSdPlot(assay(rldLiver))

pdf("plots/Liver_tRNAPhe_PCAplot.pdf", width=10, height=10)
plotPCA(vsdLiver, intgroup=c("Genotype"))
dev.off()

# NormCounts
normCountsLiver <- counts(ddsLiver, normalized=TRUE)
write.table(normCountsLiver,"tables/Liver_tRNAPhe_normCounts.txt", sep="\t", col.names=T, row.names=T, quote=F)

pdf("plots/Liver_tRNAPhe_MAplot.Phe_results.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11_Liver, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Liver, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe31_Liver, main="Phe31 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()


####################### with threshold #######################
### The effect of the Phe11 genotype - no threshold - s < 0.01 (svalues)
resPhe11_Liver.thresh <- results(ddsLiver, name="Genotype_Phe11_vs_WT", alpha=0.01)
resPhe11_Liver.thresh <- lfcShrink(ddsLiver, coef="Genotype_Phe11_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe11_Liver.thresh, svalue=TRUE)
write.table(resPhe11_Liver.thresh,"tables/Liver_Phe11_vs_WT.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11_Liver.thresh)

resPhe21_Liver.thresh <- results(ddsLiver, name="Genotype_Phe21_vs_WT", alpha=0.01)
resPhe21_Liver.thresh <- lfcShrink(ddsLiver, coef="Genotype_Phe21_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe21_Liver.thresh, svalue=TRUE)
write.table(resPhe21_Liver.thresh,"tables/Liver_Phe21_vs_WT.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21_Liver.thresh)

resPhe31_Liver.thresh <- results(ddsLiver, name="Genotype_Phe31_vs_WT", alpha=0.01)
resPhe31_Liver.thresh <- lfcShrink(ddsLiver, coef="Genotype_Phe31_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe31_Liver.thresh, svalue=TRUE)
write.table(resPhe31_Liver.thresh,"tables/Liver_Phe31_vs_WT.results.s0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe31_Liver.thresh)

pdf("plots/Liver_tRNAPhe_MAplot.Phe_results.thresh.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11_Liver.thresh, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Liver.thresh, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe31_Liver.thresh, main="Phe31 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()



#########################
### Inter-KO comparisons
#########################

#########################
# no threshold
#########################

## Phe11 vs Phe21
ddsLiver$Genotype <- relevel(ddsLiver$Genotype, ref = "Phe21")
ddsLiver <- DESeq(ddsLiver)
resultsNames(ddsLiver)

resPhe11v21_Liver <- results(ddsLiver, name="Genotype_Phe11_vs_Phe21", alpha=0.01)
resPhe11v21_Liver <- lfcShrink(ddsLiver, coef="Genotype_Phe11_vs_Phe21", type="apeglm", lfcThreshold=0, res=resPhe11v21_Liver, svalue=TRUE)
write.table(resPhe11v21_Liver,"tables/Liver_Phe11_vs_Phe21.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11v21_Liver)

## Phe11 vs Phe31
ddsLiver$Genotype <- relevel(ddsLiver$Genotype, ref = "Phe31")
ddsLiver <- DESeq(ddsLiver)
resultsNames(ddsLiver)

resPhe11v31_Liver <- results(ddsLiver, name="Genotype_Phe11_vs_Phe31", alpha=0.01)
resPhe11v31_Liver <- lfcShrink(ddsLiver, coef="Genotype_Phe11_vs_Phe31", type="apeglm", lfcThreshold=0, res=resPhe11v31_Liver, svalue=TRUE)
write.table(resPhe11v31_Liver,"tables/Liver_Phe11_vs_Phe31.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11v31_Liver)

## Phe21 vs Phe31
resPhe21v31_Liver <- results(ddsLiver, name="Genotype_Phe21_vs_Phe31", alpha=0.01)
resPhe21v31_Liver <- lfcShrink(ddsLiver, coef="Genotype_Phe21_vs_Phe31", type="apeglm", lfcThreshold=0, res=resPhe21v31_Liver, svalue=TRUE)
write.table(resPhe21v31_Liver,"tables/Liver_Phe21_vs_Phe31.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21v31_Liver)

pdf("plots/Liver_tRNAPhe_MAplot.Phe_KO_comparisons.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11v21_Liver, main="Phe11 vs Phe21", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe11vs31_Liver, main="Phe21 vs Phe31", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21vs21_Liver, main="Phe21 vs Phe31", ylim=c(-10,10), alpha=0.01)
dev.off()

#########################
# lgo2(1.5) threshold
#########################

## Phe11 vs Phe21
ddsLiver$Genotype <- relevel(ddsLiver$Genotype, ref = "Phe21")
ddsLiver <- DESeq(ddsLiver)
resultsNames(ddsLiver)

resPhe11v21_Liver.thresh <- results(ddsLiver, name="Genotype_Phe11_vs_Phe21", alpha=0.01)
resPhe11v21_Liver.thresh <- lfcShrink(ddsLiver, coef="Genotype_Phe11_vs_Phe21", type="apeglm", lfcThreshold=log2(1.5), res=resPhe11v21_Liver.thresh, svalue=TRUE)
write.table(resPhe11v21_Liver.thresh,"tables/Liver_Phe11_vs_Phe21.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11v21_Liver.thresh)

## Phe11 vs Phe31
ddsLiver$Genotype <- relevel(ddsLiver$Genotype, ref = "Phe31")
ddsLiver <- DESeq(ddsLiver)
resultsNames(ddsLiver)

resPhe11v31_Liver.thresh <- results(ddsLiver, name="Genotype_Phe11_vs_Phe31", alpha=0.01)
resPhe11v31_Liver.thresh <- lfcShrink(ddsLiver, coef="Genotype_Phe11_vs_Phe31", type="apeglm", lfcThreshold=log2(1.5), res=resPhe11v31_Liver.thresh, svalue=TRUE)
write.table(resPhe11v31_Liver.thresh,"tables/Liver_Phe11_vs_Phe31.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11v31_Liver.thresh)

## Phe21 vs Phe31
resPhe21v31_Liver.thresh <- results(ddsLiver, name="Genotype_Phe21_vs_Phe31", alpha=0.01)
resPhe21v31_Liver.thresh <- lfcShrink(ddsLiver, coef="Genotype_Phe21_vs_Phe31", type="apeglm", lfcThreshold=log2(1.5), res=resPhe21v31_Liver.thresh, svalue=TRUE)
write.table(resPhe21v31_Liver.thresh,"tables/Liver_Phe21_vs_Phe31.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21v31_Liver.thresh)

pdf("plots/Liver_tRNAPhe_MAplot.Phe_KO_comparisons.thresh.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11v21_Liver.thresh, main="Phe11 vs Phe21", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe11vs31_Liver.thresh, main="Phe21 vs Phe31", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21vs21_Liver.thresh, main="Phe21 vs Phe31", ylim=c(-10,10), alpha=0.01)
dev.off()

########################################################################################################################################################
### Brain ###
########################################################################################################################################################
samplesBrain <- subset(samples, Tissue=="Brain")
txiBrain <- tximport(files[13:24], type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")

ddsBrain <- DESeqDataSetFromTximport(txiBrain,
                                     colData = samplesBrain,
                                     design = ~ Genotype)
ddsBrain$Genotype <- relevel(ddsBrain$Genotype, ref = "WT")
ddsBrain <- DESeq(ddsBrain)
resultsNames(ddsBrain)

### The effect of the Phe11 genotype - no threshold - s < 0.01 (svalues)
resPhe11_Brain <- results(ddsBrain, name="Genotype_Phe11_vs_WT", alpha=0.01)
resPhe11_Brain <- lfcShrink(ddsBrain, coef="Genotype_Phe11_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe11_Brain, svalue=TRUE)
write.table(resPhe11_Brain,"tables/Brain_Phe11_vs_WT.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11_Brain)

resPhe21_Brain <- results(ddsBrain, name="Genotype_Phe21_vs_WT", alpha=0.01)
resPhe21_Brain <- lfcShrink(ddsBrain, coef="Genotype_Phe21_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe21_Brain, svalue=TRUE)
write.table(resPhe21_Brain,"tables/Brain_Phe21_vs_WT.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21_Brain)

resPhe31_Brain <- results(ddsBrain, name="Genotype_Phe31_vs_WT", alpha=0.01)
resPhe31_Brain <- lfcShrink(ddsBrain, coef="Genotype_Phe31_vs_WT", type="apeglm", lfcThreshold=0, res=resPhe31_Brain, svalue=TRUE)
write.table(resPhe31_Brain,"tables/Brain_Phe31_vs_WT.results.s0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe31_Brain)


## Dispersion estimates
pdf("plots/Brain_tRNAPhe_DispEstimates.pdf", width=10, height=10)
plotDispEsts(ddsBrain)
dev.off()

ntdBrain <- normTransform(ddsBrain)
rldBrain <- rlog(ddsBrain, blind=FALSE)
vsdBrain <- vst(ddsBrain, blind=FALSE)

meanSdPlot(assay(ntdBrain))
meanSdPlot(assay(vsdBrain))
meanSdPlot(assay(rldBrain))

pdf("plots/Brain_tRNAPhe_PCAplot.pdf", width=10, height=10)
plotPCA(vsdBrain, intgroup=c("Genotype"))
dev.off()

# NormCounts
normCountsBrain <- counts(ddsBrain, normalized=TRUE)
write.table(normCountsBrain,"tables/Brain_tRNAPhe_normCounts.txt", sep="\t", col.names=T, row.names=T, quote=F)

pdf("plots/Brain_tRNAPhe_MAplot.Phe11_results.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11_Brain, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Brain, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe31_Brain, main="Phe31 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()


####################### with threshold #######################
### The effect of the Phe11 genotype - no threshold - s < 0.01 (svalues)
resPhe11_Brain.thresh <- results(ddsBrain, name="Genotype_Phe11_vs_WT", alpha=0.01)
resPhe11_Brain.thresh <- lfcShrink(ddsBrain, coef="Genotype_Phe11_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe11_Brain.thresh, svalue=TRUE)
write.table(resPhe11_Brain.thresh,"tables/Brain_Phe11_vs_WT.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11_Brain.thresh)

resPhe21_Brain.thresh <- results(ddsBrain, name="Genotype_Phe21_vs_WT", alpha=0.01)
resPhe21_Brain.thresh <- lfcShrink(ddsBrain, coef="Genotype_Phe21_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe21_Brain.thresh, svalue=TRUE)
write.table(resPhe21_Brain.thresh,"tables/Brain_Phe21_vs_WT.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21_Brain.thresh)

resPhe31_Brain.thresh <- results(ddsBrain, name="Genotype_Phe31_vs_WT", alpha=0.01)
resPhe31_Brain.thresh <- lfcShrink(ddsBrain, coef="Genotype_Phe31_vs_WT", type="apeglm", lfcThreshold=log2(1.5), res=resPhe31_Brain.thresh, svalue=TRUE)
write.table(resPhe31_Brain.thresh,"tables/Brain_Phe31_vs_WT.results.s0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe31_Brain.thresh)

pdf("plots/Brain_tRNAPhe_MAplot.Phe_results.thresh.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11_Brain.thresh, main="Phe11 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21_Brain.thresh, main="Phe21 vs WT", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe31_Brain.thresh, main="Phe31 vs WT", ylim=c(-10,10), alpha=0.01)
dev.off()



#########################
### Inter-KO comparisons
#########################

#########################
# no threshold
#########################

## Phe11 vs Phe21
ddsBrain$Genotype <- relevel(ddsBrain$Genotype, ref = "Phe21")
ddsBrain <- DESeq(ddsBrain)
resultsNames(ddsBrain)

resPhe11v21_Brain <- results(ddsBrain, name="Genotype_Phe11_vs_Phe21", alpha=0.01)
resPhe11v21_Brain <- lfcShrink(ddsBrain, coef="Genotype_Phe11_vs_Phe21", type="apeglm", lfcThreshold=0, res=resPhe11v21_Brain, svalue=TRUE)
write.table(resPhe11v21_Brain,"tables/Brain_Phe11_vs_Phe21.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11v21_Brain)

## Phe11 vs Phe31
ddsBrain$Genotype <- relevel(ddsBrain$Genotype, ref = "Phe31")
ddsBrain <- DESeq(ddsBrain)
resultsNames(ddsBrain)

resPhe11v31_Brain <- results(ddsBrain, name="Genotype_Phe11_vs_Phe31", alpha=0.01)
resPhe11v31_Brain <- lfcShrink(ddsBrain, coef="Genotype_Phe11_vs_Phe31", type="apeglm", lfcThreshold=0, res=resPhe11v31_Brain, svalue=TRUE)
write.table(resPhe11v31_Brain,"tables/Brain_Phe11_vs_Phe31.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11v31_Brain)

## Phe21 vs Phe31
resPhe21v31_Brain <- results(ddsBrain, name="Genotype_Phe21_vs_Phe31", alpha=0.01)
resPhe21v31_Brain <- lfcShrink(ddsBrain, coef="Genotype_Phe21_vs_Phe31", type="apeglm", lfcThreshold=0, res=resPhe21v31_Brain, svalue=TRUE)
write.table(resPhe21v31_Brain,"tables/Brain_Phe21_vs_Phe31.results.s_0.01.no_thresh.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21v31_Brain)

pdf("plots/Brain_tRNAPhe_MAplot.Phe_KO_comparisons.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11v21_Brain, main="Phe11 vs Phe21", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe11vs31_Brain, main="Phe21 vs Phe31", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21vs21_Brain, main="Phe21 vs Phe31", ylim=c(-10,10), alpha=0.01)
dev.off()


#########################
# log2(1.5) threshold
#########################

## Phe11 vs Phe21
ddsBrain$Genotype <- relevel(ddsBrain$Genotype, ref = "Phe21")
ddsBrain <- DESeq(ddsBrain)
resultsNames(ddsBrain)

resPhe11v21_Brain.thresh <- results(ddsBrain, name="Genotype_Phe11_vs_Phe21", alpha=0.01)
resPhe11v21_Brain.thresh <- lfcShrink(ddsBrain, coef="Genotype_Phe11_vs_Phe21", type="apeglm", lfcThreshold=log2(1.5), res=resPhe11v21_Brain.thresh, svalue=TRUE)
write.table(resPhe11v21_Brain.thresh,"tables/Brain_Phe11_vs_Phe21.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11v21_Brain.thresh)

## Phe11 vs Phe31
ddsBrain$Genotype <- relevel(ddsBrain$Genotype, ref = "Phe31")
ddsBrain <- DESeq(ddsBrain)
resultsNames(ddsBrain)

resPhe11v31_Brain.thresh <- results(ddsBrain, name="Genotype_Phe11_vs_Phe31", alpha=0.01)
resPhe11v31_Brain.thresh <- lfcShrink(ddsBrain, coef="Genotype_Phe11_vs_Phe31", type="apeglm", lfcThreshold=log2(1.5), res=resPhe11v31_Brain.thresh, svalue=TRUE)
write.table(resPhe11v31_Brain.thresh,"tables/Brain_Phe11_vs_Phe31.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe11v31_Brain.thresh)

## Phe21 vs Phe31
resPhe21v31_Brain.thresh <- results(ddsBrain, name="Genotype_Phe21_vs_Phe31", alpha=0.01)
resPhe21v31_Brain.thresh <- lfcShrink(ddsBrain, coef="Genotype_Phe21_vs_Phe31", type="apeglm", lfcThreshold=log2(1.5), res=resPhe21v31_Brain.thresh, svalue=TRUE)
write.table(resPhe21v31_Brain.thresh,"tables/Brain_Phe21_vs_Phe31.results.s_0.01.lfc_log2_1.5.txt", sep="\t", col.names=T, row.names=T, quote=F)
summary(resPhe21v31_Brain.thresh)

pdf("plots/Brain_tRNAPhe_MAplot.Phe_KO_comparisons.thresh.pdf", width=15, height=10)
par(mfrow=c(1,3))
plotMA(resPhe11v21_Brain.thresh, main="Phe11 vs Phe21", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe11vs31_Brain.thresh, main="Phe21 vs Phe31", ylim=c(-10,10), alpha=0.01)
plotMA(resPhe21vs21_Brain.thresh, main="Phe21 vs Phe31", ylim=c(-10,10), alpha=0.01)
dev.off()

save.image()

####################################################################################################################################

