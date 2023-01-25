
setwd("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis")

library(DiffBind)
library(rtracklayer)

samples <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis/Phe31_WT.brain.samples.txt",header=TRUE,sep="\t")

#### brainPhe31 SAMPLE ANALYSIS - will aid in determining quality of results - PCA all vs BD - account for tissue ####
## 1. Read in brainPhe31_peaksets
brainPhe31_peaks <- dba(sampleSheet=samples)
brainPhe31_peaks

plot(brainPhe31_peaks)

## 1a. Blacklists (and greylist?)
brainPhe31_peakdata <- dba.show(brainPhe31_peaks)$Intervals
brainPhe31_peaks <- dba.blacklist(brainPhe31_peaks, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE)
brainPhe31_peaks

brainPhe31_peakdata.BL <- dba.show(brainPhe31_peaks)$Intervals
brainPhe31_peakdata - brainPhe31_peakdata.BL

## 1b. Consensus peaks
brainPhe31_consensus <- dba.peakset(brainPhe31_peaks, consensus=DBA_TISSUE, minOverlap=1)	# for replicate consensus, but has already been done by genrich
brainPhe31_consensus <- dba(brainPhe31_consensus, mask=brainPhe31_consensus$masks$"Consensus", minOverlap=1)
brainPhe31_consensus_peaks <- dba.peakset(brainPhe31_consensus, bRetrieve=TRUE)
brainPhe31_consensus_peaks

write.table(brainPhe31_consensus_peaks, "brainPhe31_v_WT_consensus_peaks.high_conf_peaks.nfr.txt", sep="\t")

## 2. Counting reads
brainPhe31_peaks <- dba.count(brainPhe31_peaks, peaks=brainPhe31_consensus_peaks, bParallel=FALSE, summits=FALSE, mapQCth=30, fragmentSize=0)
brainPhe31_peaks

brainPhe31_info <- dba.show(brainPhe31_peaks)
brainPhe31_libsizes <- cbind(LibReads=brainPhe31_info$Reads, FRiP=brainPhe31_info$FRiP, PeakReads=round(brainPhe31_info$Reads * brainPhe31_info$FRiP))
rownames(brainPhe31_libsizes) <- brainPhe31_info$ID
brainPhe31_libsizes

pdf("brainPhe31/brainPhe31_peaks.high_conf_peaks.all_sites.nfr.pdf")
plot(brainPhe31_peaks)
dev.off()

## 3. Normalisation (try all three)
brainPhe31_peaks <- dba.normalize(brainPhe31_peaks, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE) ## by default library is set to RiP
brainPhe31_norm.DESeq2 <- dba.normalize(brainPhe31_peaks, bRetrieve=TRUE, method=DBA_DESEQ2)
brainPhe31_norm.DESeq2

brainPhe31_norm.edgeR <- dba.normalize(brainPhe31_peaks, bRetrieve=TRUE, method=DBA_EDGER)
brainPhe31_norm.edgeR

brainPhe31_normlibs.DESeq2 <- cbind(FullLibSize=brainPhe31_norm.DESeq2$lib.sizes, NormFacs=brainPhe31_norm.DESeq2$norm.factors, NormLibSize=round(brainPhe31_norm.DESeq2$lib.sizes/brainPhe31_norm.DESeq2$norm.factors))
rownames(brainPhe31_normlibs.DESeq2) <- brainPhe31_info$ID
brainPhe31_normlibs.DESeq2

brainPhe31_normlibs.edgeR <- cbind(FullLibSize=brainPhe31_norm.edgeR$lib.sizes, NormFacs=brainPhe31_norm.edgeR$norm.factors, NormLibSize=round(brainPhe31_norm.edgeR$lib.sizes/brainPhe31_norm.edgeR$norm.factors))
rownames(brainPhe31_normlibs.edgeR) <- brainPhe31_info$ID
brainPhe31_normlibs.edgeR

write.table(brainPhe31_normlibs.DESeq2, "brainPhe31_normlibs.DESeq2.RiP.nfr.txt", sep="\t", quote=FALSE)
write.table(brainPhe31_normlibs.edgeR, "brainPhe31_normlibs_edgeR.RiP.nfr.txt", sep="\t", quote=FALSE)

## 4. Establishing a model design
brainPhe31_peaks <- dba.contrast(brainPhe31_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
brainPhe31_peaks

## 5. Performing the analysis
brainPhe31_peaks <- dba.analyze(brainPhe31_peaks, method=DBA_ALL_METHODS, bGreylist=FALSE)
dba.show(brainPhe31_peaks, bContrasts=TRUE)
brainPhe31_results_DESeq2 <- dba.analyze(brainPhe31_peaks, bRetrieveAnalysis=DBA_DESEQ2)
brainPhe31_results_edgeR <- dba.analyze(brainPhe31_peaks, bRetrieveAnalysis=DBA_EDGER)


pdf("brainPhe31/brainPhe31_peaks.high_conf_peaks.DESeq2.DB_peaks.nfr.pdf")
plot(brainPhe31_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe31/brainPhe31_peaks.high_conf_peaks.edgeR.DB_peaks.nfr.pdf")
plot(brainPhe31_peaks, contrast=1, method=DBA_EDGER)
dev.off()

## 6. Retrieving the DB sites
brainPhe31_peaks.DB.DESeq2 <- dba.report(brainPhe31_peaks, method=DBA_DESEQ2)
sum(brainPhe31_peaks.DB.DESeq2$Fold>0)
sum(brainPhe31_peaks.DB.DESeq2$Fold<0)
write.table(brainPhe31_peaks.DB.DESeq2, "brainPhe31_v_WT.high_conf_peaks.DESEQ2_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")

brainPhe31_peaks.DB.edgeR <- dba.report(brainPhe31_peaks, method=DBA_EDGER)
sum(brainPhe31_peaks.DB.edgeR$Fold>0)
sum(brainPhe31_peaks.DB.edgeR$Fold<0)
write.table(brainPhe31_peaks.DB.edgeR, "brainPhe31_v_WT.high_conf_peaks.EDGER_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")

brainPhe31_peaks.ALL.DESeq2 <- dba.report(brainPhe31_peaks, method=DBA_DESEQ2, th=1)
brainPhe31_peaks.ALL.edgeR <- dba.report(brainPhe31_peaks, method=DBA_EDGER, th=1)

write.table(brainPhe31_peaks.ALL.DESeq2, "brainPhe31_v_WT.high_conf_peaks.DESEQ2_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")
write.table(brainPhe31_peaks.ALL.edgeR, "brainPhe31_v_WT.high_conf_peaks.EDGER_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")

save.image("brainPhe31/brainPhe31_v_WT.high_conf_peaks.RData")

writeLines(capture.output(sessionInfo()), "sessionInfo.nfr.txt")

##############################
##			PLOTS			##
##############################

# Venn diagrams
dba.plotVenn(brainPhe31_peaks, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)

# PCA Plots
pdf("brainPhe31/brainPhe31_v_WT.high_conf_peaks.ALLSITES.PCA.nfr.pdf")
dba.plotPCA(brainPhe31_peaks, DBA_TISSUE, label=DBA_CONDITION)
dev.off()
pdf("brainPhe31/brainPhe31_v_WT.high_conf_peaks.DESEQ2.PCA.nfr.pdf")
dba.plotPCA(brainPhe31_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe31/brainPhe31_v_WT.high_conf_peaks.EDGER.PCA.nfr.pdf")
dba.plotPCA(brainPhe31_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_EDGER)
dev.off()

# MA plots
pdf("brainPhe31/brainPhe31_v_WT.high_conf_peaks.DESEQ2.MAplot.nfr.pdf")
dba.plotMA(brainPhe31_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe31/brainPhe31_v_WT.high_conf_peaks.EDGER.MAplot.nfr.pdf")
dba.plotMA(brainPhe31_peaks, contrast=1, method=DBA_EDGER)
dev.off()

# Volcano plots
dba.plotVolcano(brainPhe31_peaks, contrast=1, method=DBA_DESEQ2)

# Boxplots
sum(brainPhe31_peaks.DB$Fold<0)
sum(brainPhe31_peaks.DB$Fold>0)

brainPhe31_pvals <- dba.plotBox(brainPhe31_peaks)
brainPhe31_pvals

# Heatmaps
brainPhe31_corvals <- dba.plotHeatmap(brainPhe31_peaks)
brainPhe31_hmap <- colorRampPalette(c("red","black","green"))(n=13)
brainPhe31_readscores <- dba.plotHeatmap(brainPhe31_peaks, contrast=1, correlations=FALSE, scale="row", colScheme=brainPhe31_hmap)


### LOESS NORMALISED TEST ###
library(csaw)
loess.brainPhe31_peaks <- brainPhe31_peaks
loess.brainPhe31_peaks$config$AnalysisMethod <- DBA_DESEQ2
loess.brainPhe31_peaks <- dba.normalize(loess.brainPhe31_peaks, offsets=TRUE)
loess.brainPhe31_peaks <- dba.contrast(loess.brainPhe31_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
loess.brainPhe31_peaks <- dba.analyze(loess.brainPhe31_peaks, bGreylist=FALSE)

pdf("loess.brainPhe31_peaks.EDGER.MAplot.nfr.pdf")
dba.plotMA(loess.brainPhe31_peaks)
dev.off()
loess.brainPhe31_peaks.DB.EDGER <- dba.report(loess.brainPhe31_peaks)
sum(loess.brainPhe31_peaks.DB.EDGER$Fold > 0)
sum(loess.brainPhe31_peaks.DB.EDGER$Fold < 0)

pdf("loess.brainPhe31_peaks.DESEQ2.MAplot.nfr.pdf")
dba.plotMA(loess.brainPhe31_peaks)
dev.off()
loess.brainPhe31_peaks.DB.DESEQ2 <- dba.report(loess.brainPhe31_peaks)
sum(loess.brainPhe31_peaks.DB.DESEQ2$Fold > 0)
sum(loess.brainPhe31_peaks.DB.DESEQ2$Fold < 0)

write.table(loess.brainPhe31_peaks.DB.EDGER, "brainPhe31_v_WT.high_conf_peaks.loess_EDGER.DB_results.nfr.txt", quote=FALSE, sep="\t")
write.table(loess.brainPhe31_peaks.DB.DESEQ2, "brainPhe31_v_WT.high_conf_peaks.loess_DESEQ2.DB_results.nfr.txt", quote=FALSE, sep="\t")


writeLines(capture.output(brainPhe31_peaks), "brainPhe31_peaks.nfr.txt")


