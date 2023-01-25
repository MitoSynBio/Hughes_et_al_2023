
setwd("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc")

library(DiffBind)
library(rtracklayer)

samples <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/diffbind_samples_filt.redux.txt",header=TRUE,sep="\t")

liverPhe31Samples <- subset(samples, Tissue=="Liver" & Condition=="WT" | Tissue=="Liver" & Condition=="Phe31")

#### liverPhe31 SAMPLE ANALYSIS - will aid in determining quality of results - PCA all vs BD - account for tissue ####
## 1. Read in liverPhe31_peaksets
liverPhe31_peaks <- dba(sampleSheet=liverPhe31Samples)
liverPhe31_peaks

plot(liverPhe31_peaks)

## 1a. Blacklists (and greylist?)
liverPhe31_peakdata <- dba.show(liverPhe31_peaks)$Intervals
liverPhe31_peaks <- dba.blacklist(liverPhe31_peaks, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE)
liverPhe31_peaks

liverPhe31_peakdata.BL <- dba.show(liverPhe31_peaks)$Intervals
liverPhe31_peakdata - liverPhe31_peakdata.BL

## 1b. Consensus peaks
liverPhe31_consensus <- dba.peakset(liverPhe31_peaks, consensus=DBA_TISSUE, minOverlap=1)	# for replicate consensus, but has already been done by genrich
liverPhe31_consensus <- dba(liverPhe31_consensus, mask=liverPhe31_consensus$masks$"Consensus", minOverlap=1)
liverPhe31_consensus_peaks <- dba.peakset(liverPhe31_consensus, bRetrieve=TRUE)
liverPhe31_consensus_peaks

write.table(liverPhe31_consensus_peaks, "redux.liverPhe31/liverPhe31_v_WT_consensus_peaks.high_conf_peaks.txt", sep="\t")

## 2. Counting reads
liverPhe31_peaks <- dba.count(liverPhe31_peaks, peaks=liverPhe31_consensus_peaks, bParallel=FALSE, summits=FALSE, mapQCth=30, fragmentSize=0)
liverPhe31_peaks

liverPhe31_info <- dba.show(liverPhe31_peaks)
liverPhe31_libsizes <- cbind(LibReads=liverPhe31_info$Reads, FRiP=liverPhe31_info$FRiP, PeakReads=round(liverPhe31_info$Reads * liverPhe31_info$FRiP))
rownames(liverPhe31_libsizes) <- liverPhe31_info$ID
liverPhe31_libsizes

pdf("redux.liverPhe31/liverPhe31_peaks.high_conf_peaks.all_sites.pdf")
plot(liverPhe31_peaks)
dev.off()

## 3. Normalisation (try all three)
liverPhe31_peaks <- dba.normalize(liverPhe31_peaks, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE) ## by default library is set to RiP
liverPhe31_norm.DESeq2 <- dba.normalize(liverPhe31_peaks, bRetrieve=TRUE, method=DBA_DESEQ2)
liverPhe31_norm.DESeq2

liverPhe31_norm.edgeR <- dba.normalize(liverPhe31_peaks, bRetrieve=TRUE, method=DBA_EDGER)
liverPhe31_norm.edgeR

liverPhe31_normlibs.DESeq2 <- cbind(FullLibSize=liverPhe31_norm.DESeq2$lib.sizes, NormFacs=liverPhe31_norm.DESeq2$norm.factors, NormLibSize=round(liverPhe31_norm.DESeq2$lib.sizes/liverPhe31_norm.DESeq2$norm.factors))
rownames(liverPhe31_normlibs.DESeq2) <- liverPhe31_info$ID
liverPhe31_normlibs.DESeq2

liverPhe31_normlibs.edgeR <- cbind(FullLibSize=liverPhe31_norm.edgeR$lib.sizes, NormFacs=liverPhe31_norm.edgeR$norm.factors, NormLibSize=round(liverPhe31_norm.edgeR$lib.sizes/liverPhe31_norm.edgeR$norm.factors))
rownames(liverPhe31_normlibs.edgeR) <- liverPhe31_info$ID
liverPhe31_normlibs.edgeR

write.table(liverPhe31_normlibs.DESeq2, "redux.liverPhe31/liverPhe31_normlibs.DESeq2.RiP.txt", sep="\t", quote=FALSE)
write.table(liverPhe31_normlibs.edgeR, "redux.liverPhe31/liverPhe31_normlibs_edgeR.RiP.txt", sep="\t", quote=FALSE)

## 4. Establishing a model design
liverPhe31_peaks <- dba.contrast(liverPhe31_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
liverPhe31_peaks

## 5. Performing the analysis
liverPhe31_peaks <- dba.analyze(liverPhe31_peaks, method=DBA_ALL_METHODS, bGreylist=FALSE)
dba.show(liverPhe31_peaks, bContrasts=TRUE)
liverPhe31_results_DESeq2 <- dba.analyze(liverPhe31_peaks, bRetrieveAnalysis=DBA_DESEQ2)
liverPhe31_results_edgeR <- dba.analyze(liverPhe31_peaks, bRetrieveAnalysis=DBA_EDGER)


pdf("redux.liverPhe31/liverPhe31_peaks.high_conf_peaks.DESeq2.DB_peaks.pdf")
plot(liverPhe31_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("redux.liverPhe31/liverPhe31_peaks.high_conf_peaks.edgeR.DB_peaks.pdf")
plot(liverPhe31_peaks, contrast=1, method=DBA_EDGER)
dev.off()

## 6. Retrieving the DB sites
liverPhe31_peaks.DB.DESeq2 <- dba.report(liverPhe31_peaks, method=DBA_DESEQ2)
sum(liverPhe31_peaks.DB.DESeq2$Fold>0)
sum(liverPhe31_peaks.DB.DESeq2$Fold<0)
write.table(liverPhe31_peaks.DB.DESeq2, "redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.DESEQ2_RiP.DB_results.txt", quote=FALSE, sep="\t")

liverPhe31_peaks.DB.edgeR <- dba.report(liverPhe31_peaks, method=DBA_EDGER)
sum(liverPhe31_peaks.DB.edgeR$Fold>0)
sum(liverPhe31_peaks.DB.edgeR$Fold<0)
write.table(liverPhe31_peaks.DB.edgeR, "redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.EDGER_RiP.DB_results.txt", quote=FALSE, sep="\t")

liverPhe31_peaks.ALL.DESeq2 <- dba.report(liverPhe31_peaks, method=DBA_DESEQ2, th=1)
liverPhe31_peaks.ALL.edgeR <- dba.report(liverPhe31_peaks, method=DBA_EDGER, th=1)

write.table(liverPhe31_peaks.ALL.DESeq2, "redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.DESEQ2_RiP.ALL_results.txt", quote=FALSE, sep="\t")
write.table(liverPhe31_peaks.ALL.edgeR, "redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.EDGER_RiP.ALL_results.txt", quote=FALSE, sep="\t")

save.image("redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.RData")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

##############################
##			PLOTS			##
##############################

# Venn diagrams
dba.plotVenn(liverPhe31_peaks, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)

# PCA Plots
pdf("redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.ALLSITES.PCA.pdf")
dba.plotPCA(liverPhe31_peaks, DBA_TISSUE, label=DBA_CONDITION)
dev.off()
pdf("redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.DESEQ2.PCA.pdf")
dba.plotPCA(liverPhe31_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.EDGER.PCA.pdf")
dba.plotPCA(liverPhe31_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_EDGER)
dev.off()

# MA plots
pdf("redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.DESEQ2.MAplot.pdf")
dba.plotMA(liverPhe31_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.EDGER.MAplot.pdf")
dba.plotMA(liverPhe31_peaks, contrast=1, method=DBA_EDGER)
dev.off()

# Volcano plots
dba.plotVolcano(liverPhe31_peaks, contrast=1, method=DBA_DESEQ2)

# Boxplots
sum(liverPhe31_peaks.DB$Fold<0)
sum(liverPhe31_peaks.DB$Fold>0)

liverPhe31_pvals <- dba.plotBox(liverPhe31_peaks)
liverPhe31_pvals

# Heatmaps
liverPhe31_corvals <- dba.plotHeatmap(liverPhe31_peaks)
liverPhe31_hmap <- colorRampPalette(c("red","black","green"))(n=13)
liverPhe31_readscores <- dba.plotHeatmap(liverPhe31_peaks, contrast=1, correlations=FALSE, scale="row", colScheme=liverPhe31_hmap)


### LOESS NORMALISED TEST ###
library(csaw)
loess.liverPhe31_peaks <- liverPhe31_peaks
loess.liverPhe31_peaks$config$AnalysisMethod <- DBA_DESEQ2
loess.liverPhe31_peaks <- dba.normalize(loess.liverPhe31_peaks, offsets=TRUE)
loess.liverPhe31_peaks <- dba.contrast(loess.liverPhe31_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
loess.liverPhe31_peaks <- dba.analyze(loess.liverPhe31_peaks, bGreylist=FALSE)

pdf("redux.liverPhe31/loess.liverPhe31_peaks.EDGER.MAplot.pdf")
dba.plotMA(loess.liverPhe31_peaks)
dev.off()
loess.liverPhe31_peaks.DB.EDGER <- dba.report(loess.liverPhe31_peaks)
sum(loess.liverPhe31_peaks.DB.EDGER$Fold > 0)
sum(loess.liverPhe31_peaks.DB.EDGER$Fold < 0)

pdf("redux.liverPhe31/loess.liverPhe31_peaks.DESEQ2.MAplot.pdf")
dba.plotMA(loess.liverPhe31_peaks)
dev.off()
loess.liverPhe31_peaks.DB.DESEQ2 <- dba.report(loess.liverPhe31_peaks)
sum(loess.liverPhe31_peaks.DB.DESEQ2$Fold > 0)
sum(loess.liverPhe31_peaks.DB.DESEQ2$Fold < 0)

write.table(loess.liverPhe31_peaks.DB.EDGER, "redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.loess_EDGER.DB_results.txt", quote=FALSE, sep="\t")
write.table(loess.liverPhe31_peaks.DB.DESEQ2, "redux.liverPhe31/liverPhe31_v_WT.high_conf_peaks.loess_DESEQ2.DB_results.txt", quote=FALSE, sep="\t")


writeLines(capture.output(liverPhe31_peaks), "redux.liverPhe31/liverPhe31_peaks.txt")


