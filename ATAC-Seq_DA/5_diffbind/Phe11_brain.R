
setwd("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis")

library(DiffBind)
library(rtracklayer)

samples <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis/Phe11_WT.brain.samples.txt",header=TRUE,sep="\t")

#### brainPhe11 SAMPLE ANALYSIS - will aid in determining quality of results - PCA all vs BD - account for tissue ####
## 1. Read in brainPhe11_peaksets
brainPhe11_peaks <- dba(sampleSheet=samples)
brainPhe11_peaks

plot(brainPhe11_peaks)

## 1a. Blacklists (and greylist?)
brainPhe11_peakdata <- dba.show(brainPhe11_peaks)$Intervals
brainPhe11_peaks <- dba.blacklist(brainPhe11_peaks, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE)
brainPhe11_peaks

brainPhe11_peakdata.BL <- dba.show(brainPhe11_peaks)$Intervals
brainPhe11_peakdata - brainPhe11_peakdata.BL

## 1b. Consensus peaks
brainPhe11_consensus <- dba.peakset(brainPhe11_peaks, consensus=DBA_TISSUE, minOverlap=1)	# for replicate consensus, but has already been done by genrich
brainPhe11_consensus <- dba(brainPhe11_consensus, mask=brainPhe11_consensus$masks$"Consensus", minOverlap=1)
brainPhe11_consensus_peaks <- dba.peakset(brainPhe11_consensus, bRetrieve=TRUE)
brainPhe11_consensus_peaks

write.table(brainPhe11_consensus_peaks, "brainPhe11/brainPhe11_v_WT_consensus_peaks.high_conf_peaks.nfr.txt", sep="\t")

## 2. Counting reads
brainPhe11_peaks <- dba.count(brainPhe11_peaks, peaks=brainPhe11_consensus_peaks, bParallel=FALSE, summits=FALSE, mapQCth=30, fragmentSize=0)
brainPhe11_peaks

brainPhe11_info <- dba.show(brainPhe11_peaks)
brainPhe11_libsizes <- cbind(LibReads=brainPhe11_info$Reads, FRiP=brainPhe11_info$FRiP, PeakReads=round(brainPhe11_info$Reads * brainPhe11_info$FRiP))
rownames(brainPhe11_libsizes) <- brainPhe11_info$ID
brainPhe11_libsizes

pdf("brainPhe11/brainphe11_peaks.high_conf_peaks.all_sites.nfr.pdf")
plot(brainPhe11_peaks)
dev.off()

## 3. Normalisation (try all three)
brainPhe11_peaks <- dba.normalize(brainPhe11_peaks, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE) ## by default library is set to RiP
brainPhe11_norm.DESeq2 <- dba.normalize(brainPhe11_peaks, bRetrieve=TRUE, method=DBA_DESEQ2)
brainPhe11_norm.DESeq2

brainPhe11_norm.edgeR <- dba.normalize(brainPhe11_peaks, bRetrieve=TRUE, method=DBA_EDGER)
brainPhe11_norm.edgeR

brainPhe11_normlibs.DESeq2 <- cbind(FullLibSize=brainPhe11_norm.DESeq2$lib.sizes, NormFacs=brainPhe11_norm.DESeq2$norm.factors, NormLibSize=round(brainPhe11_norm.DESeq2$lib.sizes/brainPhe11_norm.DESeq2$norm.factors))
rownames(brainPhe11_normlibs.DESeq2) <- brainPhe11_info$ID
brainPhe11_normlibs.DESeq2

brainPhe11_normlibs.edgeR <- cbind(FullLibSize=brainPhe11_norm.edgeR$lib.sizes, NormFacs=brainPhe11_norm.edgeR$norm.factors, NormLibSize=round(brainPhe11_norm.edgeR$lib.sizes/brainPhe11_norm.edgeR$norm.factors))
rownames(brainPhe11_normlibs.edgeR) <- brainPhe11_info$ID
brainPhe11_normlibs.edgeR

write.table(brainPhe11_normlibs.DESeq2, "brainPhe11/brainPhe11_normlibs.DESeq2.RiP.nfr.txt", sep="\t", quote=FALSE)
write.table(brainPhe11_normlibs.edgeR, "brainPhe11/brainPhe11_normlibs_edgeR.RiP.nfr.txt", sep="\t", quote=FALSE)

## 4. Establishing a model design
brainPhe11_peaks <- dba.contrast(brainPhe11_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
brainPhe11_peaks

## 5. Performing the analysis
brainPhe11_peaks <- dba.analyze(brainPhe11_peaks, method=DBA_ALL_METHODS, bGreylist=FALSE)
dba.show(brainPhe11_peaks, bContrasts=TRUE)
brainPhe11_results_DESeq2 <- dba.analyze(brainPhe11_peaks, bRetrieveAnalysis=DBA_DESEQ2)
brainPhe11_results_edgeR <- dba.analyze(brainPhe11_peaks, bRetrieveAnalysis=DBA_EDGER)


pdf("brainPhe11/brainPhe11_peaks.high_conf_peaks.DESeq2.DB_peaks.nfr.pdf")
plot(brainPhe11_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe11/brainPhe11_peaks.high_conf_peaks.edgeR.DB_peaks.nfr.pdf")
plot(brainPhe11_peaks, contrast=1, method=DBA_EDGER)
dev.off()

## 6. Retrieving the DB sites
brainPhe11_peaks.DB.DESeq2 <- dba.report(brainPhe11_peaks, method=DBA_DESEQ2)
sum(brainPhe11_peaks.DB.DESeq2$Fold>0)
sum(brainPhe11_peaks.DB.DESeq2$Fold<0)
write.table(brainPhe11_peaks.DB.DESeq2, "brainPhe11/brainPhe11_v_WT.high_conf_peaks.DESEQ2_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")

brainPhe11_peaks.DB.edgeR <- dba.report(brainPhe11_peaks, method=DBA_EDGER)
sum(brainPhe11_peaks.DB.edgeR$Fold>0)
sum(brainPhe11_peaks.DB.edgeR$Fold<0)
write.table(brainPhe11_peaks.DB.edgeR, "brainPhe11/brainPhe11_v_WT.high_conf_peaks.EDGER_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")

brainPhe11_peaks.ALL.DESeq2 <- dba.report(brainPhe11_peaks, method=DBA_DESEQ2, th=1)
brainPhe11_peaks.ALL.edgeR <- dba.report(brainPhe11_peaks, method=DBA_EDGER, th=1)

write.table(brainPhe11_peaks.ALL.DESeq2, "brainPhe11/brainPhe11_v_WT.high_conf_peaks.DESEQ2_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")
write.table(brainPhe11_peaks.ALL.edgeR, "brainPhe11/brainPhe11_v_WT.high_conf_peaks.EDGER_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")

save.image("brainPhe11/brainPhe11_v_WT.high_conf_peaks.RData")

##############################
##			PLOTS			##
##############################

# Venn diagrams
dba.plotVenn(brainPhe11_peaks, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)

# PCA Plots
pdf("brainPhe11/brainPhe11_v_WT.high_conf_peaks.ALLSITES.PCA.nfr.pdf")
dba.plotPCA(brainPhe11_peaks, DBA_TISSUE, label=DBA_CONDITION)
dev.off()
pdf("brainPhe11/brainPhe11_v_WT.high_conf_peaks.DESEQ2.PCA.nfr.pdf")
dba.plotPCA(brainPhe11_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe11/brainPhe11_v_WT.high_conf_peaks.EDGER.PCA.nfr.pdf")
dba.plotPCA(brainPhe11_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_EDGER)
dev.off()

# MA plots
pdf("brainPhe11/brainPhe11_v_WT.high_conf_peaks.DESEQ2.MAplot.nfr.pdf")
dba.plotMA(brainPhe11_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("brainPhe11/brainPhe11_v_WT.high_conf_peaks.EDGER.MAplot.nfr.pdf")
dba.plotMA(brainPhe11_peaks, contrast=1, method=DBA_EDGER)
dev.off()

# Volcano plots
dba.plotVolcano(brainPhe11_peaks, contrast=1, method=DBA_DESEQ2)
dba.plotVolcano(brainPhe11_peaks, contrast=1, method=DBA_EDGER)

# Boxplots
sum(brainPhe11_peaks.DB$Fold<0)
sum(brainPhe11_peaks.DB$Fold>0)

brainPhe11_pvals <- dba.plotBox(brainPhe11_peaks)
brainPhe11_pvals

# Heatmaps
brainPhe11_corvals <- dba.plotHeatmap(brainPhe11_peaks)
brainPhe11_hmap <- colorRampPalette(c("red","black","green"))(n=13)
brainPhe11_readscores <- dba.plotHeatmap(brainPhe11_peaks, contrast=1, correlations=FALSE, scale="row", colScheme=brainPhe11_hmap)

