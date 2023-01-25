
setwd("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis")

library(DiffBind)
library(rtracklayer)

samples <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis/Phe21_WT.brain.samples.txt",header=TRUE,sep="\t")

#### brainPhe21 SAMPLE ANALYSIS - will aid in determining quality of results - PCA all vs BD - account for tissue ####
## 1. Read in brainPhe21_peaksets
brainPhe21_peaks <- dba(sampleSheet=samples)
brainPhe21_peaks

plot(brainPhe21_peaks)

## 1a. Blacklists (and greylist?)
brainPhe21_peakdata <- dba.show(brainPhe21_peaks)$Intervals
brainPhe21_peaks <- dba.blacklist(brainPhe21_peaks, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE)
brainPhe21_peaks

brainPhe21_peakdata.BL <- dba.show(brainPhe21_peaks)$Intervals
brainPhe21_peakdata - brainPhe21_peakdata.BL

## 1b. Consensus peaks
brainPhe21_consensus <- dba.peakset(brainPhe21_peaks, consensus=DBA_TISSUE, minOverlap=1)	# for replicate consensus, but has already been done by genrich
brainPhe21_consensus <- dba(brainPhe21_consensus, mask=brainPhe21_consensus$masks$"Consensus", minOverlap=1)
brainPhe21_consensus_peaks <- dba.peakset(brainPhe21_consensus, bRetrieve=TRUE)
brainPhe21_consensus_peaks

write.table(brainPhe21_consensus_peaks, "brainPhe21_v_WT_consensus_peaks.high_conf_peaks.nfr.txt", sep="\t")

## 2. Counting reads
brainPhe21_peaks <- dba.count(brainPhe21_peaks, peaks=brainPhe21_consensus_peaks, bParallel=FALSE, summits=FALSE, mapQCth=30, fragmentSize=0)
brainPhe21_peaks

brainPhe21_info <- dba.show(brainPhe21_peaks)
brainPhe21_libsizes <- cbind(LibReads=brainPhe21_info$Reads, FRiP=brainPhe21_info$FRiP, PeakReads=round(brainPhe21_info$Reads * brainPhe21_info$FRiP))
rownames(brainPhe21_libsizes) <- brainPhe21_info$ID
brainPhe21_libsizes

pdf("brainPhe21/brainPhe21_peaks.high_conf_peaks.all_sites.nfr.pdf")
plot(brainPhe21_peaks)
dev.off()

## 3. Normalisation (try all three)
brainPhe21_peaks <- dba.normalize(brainPhe21_peaks, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE) ## by default library is set to RiP
brainPhe21_norm.DESeq2 <- dba.normalize(brainPhe21_peaks, bRetrieve=TRUE, method=DBA_DESEQ2)
brainPhe21_norm.DESeq2

brainPhe21_norm.edgeR <- dba.normalize(brainPhe21_peaks, bRetrieve=TRUE, method=DBA_EDGER)
brainPhe21_norm.edgeR

brainPhe21_normlibs.DESeq2 <- cbind(FullLibSize=brainPhe21_norm.DESeq2$lib.sizes, NormFacs=brainPhe21_norm.DESeq2$norm.factors, NormLibSize=round(brainPhe21_norm.DESeq2$lib.sizes/brainPhe21_norm.DESeq2$norm.factors))
rownames(brainPhe21_normlibs.DESeq2) <- brainPhe21_info$ID
brainPhe21_normlibs.DESeq2

brainPhe21_normlibs.edgeR <- cbind(FullLibSize=brainPhe21_norm.edgeR$lib.sizes, NormFacs=brainPhe21_norm.edgeR$norm.factors, NormLibSize=round(brainPhe21_norm.edgeR$lib.sizes/brainPhe21_norm.edgeR$norm.factors))
rownames(brainPhe21_normlibs.edgeR) <- brainPhe21_info$ID
brainPhe21_normlibs.edgeR

write.table(brainPhe21_normlibs.DESeq2, "brainPhe21_normlibs.DESeq2.RiP.nfr.txt", sep="\t", quote=FALSE)
write.table(brainPhe21_normlibs.edgeR, "brainPhe21_normlibs_edgeR.RiP.nfr.txt", sep="\t", quote=FALSE)

## 4. Establishing a model design
brainPhe21_peaks <- dba.contrast(brainPhe21_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
brainPhe21_peaks

## 5. Performing the analysis
brainPhe21_peaks <- dba.analyze(brainPhe21_peaks, method=DBA_ALL_METHODS, bGreylist=FALSE)
dba.show(brainPhe21_peaks, bContrasts=TRUE)
brainPhe21_results_DESeq2 <- dba.analyze(brainPhe21_peaks, bRetrieveAnalysis=DBA_DESEQ2)
brainPhe21_results_edgeR <- dba.analyze(brainPhe21_peaks, bRetrieveAnalysis=DBA_EDGER)


pdf("brainPhe21/brainPhe21_peaks.high_conf_peaks.DESeq2.DB_peaks.nfr.pdf")
plot(brainPhe21_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe21/brainPhe21_peaks.high_conf_peaks.edgeR.DB_peaks.nfr.pdf")
plot(brainPhe21_peaks, contrast=1, method=DBA_EDGER)
dev.off()

## 6. Retrieving the DB sites
brainPhe21_peaks.DB.DESeq2 <- dba.report(brainPhe21_peaks, method=DBA_DESEQ2)
sum(brainPhe21_peaks.DB.DESeq2$Fold>0)
sum(brainPhe21_peaks.DB.DESeq2$Fold<0)
write.table(brainPhe21_peaks.DB.DESeq2, "brainPhe21_v_WT.high_conf_peaks.DESEQ2_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")

brainPhe21_peaks.DB.edgeR <- dba.report(brainPhe21_peaks, method=DBA_EDGER)
sum(brainPhe21_peaks.DB.edgeR$Fold>0)
sum(brainPhe21_peaks.DB.edgeR$Fold<0)
write.table(brainPhe21_peaks.DB.edgeR, "brainPhe21_v_WT.high_conf_peaks.EDGER_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")

brainPhe21_peaks.ALL.DESeq2 <- dba.report(brainPhe21_peaks, method=DBA_DESEQ2, th=1)
brainPhe21_peaks.ALL.edgeR <- dba.report(brainPhe21_peaks, method=DBA_EDGER, th=1)

write.table(brainPhe21_peaks.ALL.DESeq2, "brainPhe21_v_WT.high_conf_peaks.DESEQ2_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")
write.table(brainPhe21_peaks.ALL.edgeR, "brainPhe21_v_WT.high_conf_peaks.EDGER_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")

save.image("brainPhe21/brainPhe21_v_WT.high_conf_peaks.RData")

writeLines(capture.output(sessionInfo()), "sessionInfo.nfr.txt")

##############################
##			PLOTS			##
##############################

# Venn diagrams
dba.plotVenn(brainPhe21_peaks, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)

# PCA Plots
pdf("brainPhe21/brainPhe21_v_WT.high_conf_peaks.ALLSITES.PCA.nfr.pdf")
dba.plotPCA(brainPhe21_peaks, DBA_TISSUE, label=DBA_CONDITION)
dev.off()
pdf("brainPhe21/brainPhe21_v_WT.high_conf_peaks.DESEQ2.PCA.nfr.pdf")
dba.plotPCA(brainPhe21_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe21/brainPhe21_v_WT.high_conf_peaks.EDGER.PCA.nfr.pdf")
dba.plotPCA(brainPhe21_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_EDGER)
dev.off()

#dba.plotPCA(brainPhe21_peaks, contrast=1, label=DBA_TISSUE)
#dba.plotPCA(brainPhe21_peaks, attributes=c(DBA_TISSUE,DBA_CONDITION), label=DBA_REPLICATE)
#dba.plotPCA(brainPhe21_peaks, contrast=1, method=DBA_DESEQ2)

# MA plots
pdf("brainPhe21/brainPhe21_v_WT.high_conf_peaks.DESEQ2.MAplot.nfr.pdf")
dba.plotMA(brainPhe21_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe21/brainPhe21_v_WT.high_conf_peaks.EDGER.MAplot.nfr.pdf")
dba.plotMA(brainPhe21_peaks, contrast=1, method=DBA_EDGER)
dev.off()

# Volcano plots
dba.plotVolcano(brainPhe21_peaks, contrast=1, method=DBA_DESEQ2)
dba.plotVolcano(brainPhe21_peaks, contrast=1, method=DBA_EDGER)

# Boxplots
sum(brainPhe21_peaks.DB$Fold<0)
sum(brainPhe21_peaks.DB$Fold>0)

brainPhe21_pvals <- dba.plotBox(brainPhe21_peaks)
brainPhe21_pvals

# Heatmaps
brainPhe21_corvals <- dba.plotHeatmap(brainPhe21_peaks)
brainPhe21_hmap <- colorRampPalette(c("red","black","green"))(n=13)
brainPhe21_readscores <- dba.plotHeatmap(brainPhe21_peaks, contrast=1, correlations=FALSE, scale="row", colScheme=brainPhe21_hmap)

writeLines(capture.output(brainPhe21_peaks), "brainPhe21_peaks.nfr.txt")

