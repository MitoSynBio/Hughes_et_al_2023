
setwd("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis")

library(DiffBind)
library(rtracklayer)

samples <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis/Phe21_WT.liver.samples.txt",header=TRUE,sep="\t")

#### liverPhe21 SAMPLE ANALYSIS - will aid in determining quality of results - PCA all vs BD - account for tissue ####
## 1. Read in liverPhe21_peaksets
liverPhe21_peaks <- dba(sampleSheet=samples)
liverPhe21_peaks

plot(liverPhe21_peaks)

## 1a. Blacklists (and greylist?)
liverPhe21_peakdata <- dba.show(liverPhe21_peaks)$Intervals
liverPhe21_peaks <- dba.blacklist(liverPhe21_peaks, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE)
liverPhe21_peaks

liverPhe21_peakdata.BL <- dba.show(liverPhe21_peaks)$Intervals
liverPhe21_peakdata - liverPhe21_peakdata.BL

## 1b. Consensus peaks
liverPhe21_consensus <- dba.peakset(liverPhe21_peaks, consensus=DBA_TISSUE, minOverlap=1)	# for replicate consensus, but has already been done by genrich
liverPhe21_consensus <- dba(liverPhe21_consensus, mask=liverPhe21_consensus$masks$"Consensus", minOverlap=1)
liverPhe21_consensus_peaks <- dba.peakset(liverPhe21_consensus, bRetrieve=TRUE)
liverPhe21_consensus_peaks

write.table(liverPhe21_consensus_peaks, "liverPhe21_v_WT_consensus_peaks.high_conf_peaks.nfr.txt", sep="\t")

## 2. Counting reads
liverPhe21_peaks <- dba.count(liverPhe21_peaks, peaks=liverPhe21_consensus_peaks, bParallel=FALSE, summits=FALSE, mapQCth=30, fragmentSize=0)
liverPhe21_peaks

liverPhe21_info <- dba.show(liverPhe21_peaks)
liverPhe21_libsizes <- cbind(LibReads=liverPhe21_info$Reads, FRiP=liverPhe21_info$FRiP, PeakReads=round(liverPhe21_info$Reads * liverPhe21_info$FRiP))
rownames(liverPhe21_libsizes) <- liverPhe21_info$ID
liverPhe21_libsizes

pdf("liverPhe21/liverPhe21_peaks.high_conf_peaks.all_sites.nfr.pdf")
plot(liverPhe21_peaks)
dev.off()

## 3. Normalisation (try all three)
liverPhe21_peaks <- dba.normalize(liverPhe21_peaks, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE) ## by default library is set to RiP
liverPhe21_norm.DESeq2 <- dba.normalize(liverPhe21_peaks, bRetrieve=TRUE, method=DBA_DESEQ2)
liverPhe21_norm.DESeq2

liverPhe21_norm.edgeR <- dba.normalize(liverPhe21_peaks, bRetrieve=TRUE, method=DBA_EDGER)
liverPhe21_norm.edgeR

liverPhe21_normlibs.DESeq2 <- cbind(FullLibSize=liverPhe21_norm.DESeq2$lib.sizes, NormFacs=liverPhe21_norm.DESeq2$norm.factors, NormLibSize=round(liverPhe21_norm.DESeq2$lib.sizes/liverPhe21_norm.DESeq2$norm.factors))
rownames(liverPhe21_normlibs.DESeq2) <- liverPhe21_info$ID
liverPhe21_normlibs.DESeq2

liverPhe21_normlibs.edgeR <- cbind(FullLibSize=liverPhe21_norm.edgeR$lib.sizes, NormFacs=liverPhe21_norm.edgeR$norm.factors, NormLibSize=round(liverPhe21_norm.edgeR$lib.sizes/liverPhe21_norm.edgeR$norm.factors))
rownames(liverPhe21_normlibs.edgeR) <- liverPhe21_info$ID
liverPhe21_normlibs.edgeR

write.table(liverPhe21_normlibs.DESeq2, "liverPhe21_normlibs.DESeq2.RiP.nfr.txt", sep="\t", quote=FALSE)
write.table(liverPhe21_normlibs.edgeR, "liverPhe21_normlibs_edgeR.RiP.nfr.txt", sep="\t", quote=FALSE)

## 4. Establishing a model design
liverPhe21_peaks <- dba.contrast(liverPhe21_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
liverPhe21_peaks

## 5. Performing the analysis
liverPhe21_peaks <- dba.analyze(liverPhe21_peaks, method=DBA_ALL_METHODS, bGreylist=FALSE)
dba.show(liverPhe21_peaks, bContrasts=TRUE)
liverPhe21_results_DESeq2 <- dba.analyze(liverPhe21_peaks, bRetrieveAnalysis=DBA_DESEQ2)
liverPhe21_results_edgeR <- dba.analyze(liverPhe21_peaks, bRetrieveAnalysis=DBA_EDGER)


pdf("liverPhe21/liverPhe21_peaks.high_conf_peaks.DESeq2.DB_peaks.nfr.pdf")
plot(liverPhe21_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("liverPhe21/liverPhe21_peaks.high_conf_peaks.edgeR.DB_peaks.nfr.pdf")
plot(liverPhe21_peaks, contrast=1, method=DBA_EDGER)
dev.off()

## 6. Retrieving the DB sites
liverPhe21_peaks.DB.DESeq2 <- dba.report(liverPhe21_peaks, method=DBA_DESEQ2)
sum(liverPhe21_peaks.DB.DESeq2$Fold>0)
sum(liverPhe21_peaks.DB.DESeq2$Fold<0)
write.table(liverPhe21_peaks.DB.DESeq2, "liverPhe21_v_WT.high_conf_peaks.DESEQ2_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")

liverPhe21_peaks.DB.edgeR <- dba.report(liverPhe21_peaks, method=DBA_EDGER)
sum(liverPhe21_peaks.DB.edgeR$Fold>0)
sum(liverPhe21_peaks.DB.edgeR$Fold<0)
write.table(liverPhe21_peaks.DB.edgeR, "liverPhe21_v_WT.high_conf_peaks.EDGER_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")

liverPhe21_peaks.ALL.DESeq2 <- dba.report(liverPhe21_peaks, method=DBA_DESEQ2, th=1)
liverPhe21_peaks.ALL.edgeR <- dba.report(liverPhe21_peaks, method=DBA_EDGER, th=1)
write.table(liverPhe21_peaks.ALL.DESeq2, "liverPhe21_v_WT.high_conf_peaks.DESEQ2_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")
write.table(liverPhe21_peaks.ALL.edgeR, "liverPhe21_v_WT.high_conf_peaks.EDGER_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")

save.image("liverPhe21/liverPhe21_v_WT.high_conf_peaks.RData")

writeLines(capture.output(sessionInfo()), "sessionInfo.nfr.txt")

##############################
##			PLOTS			##
##############################

# Venn diagrams
dba.plotVenn(liverPhe21_peaks, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)

# PCA Plots

# PCA Plots
pdf("liverPhe21/liverPhe21_v_WT.high_conf_peaks.ALLSITES.PCA.nfr.pdf")
dba.plotPCA(liverPhe21_peaks, DBA_TISSUE, label=DBA_CONDITION)
dev.off()
pdf("liverPhe21/liverPhe21_v_WT.high_conf_peaks.DESEQ2.PCA.nfr.pdf")
dba.plotPCA(liverPhe21_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("liverPhe21/liverPhe21_v_WT.high_conf_peaks.EDGER.PCA.nfr.pdf")
dba.plotPCA(liverPhe21_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_EDGER)
dev.off()

#dba.plotPCA(liverPhe21_peaks, contrast=1, label=DBA_TISSUE)
#dba.plotPCA(liverPhe21_peaks, attributes=c(DBA_TISSUE,DBA_CONDITION), label=DBA_REPLICATE)
#dba.plotPCA(liverPhe21_peaks, contrast=1, method=DBA_DESEQ2)

# MA plots
pdf("liverPhe21/liverPhe21_v_WT.high_conf_peaks.DESEQ2.MAplot.nfr.pdf")
dba.plotMA(liverPhe21_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("liverPhe21/liverPhe21_v_WT.high_conf_peaks.EDGER.MAplot.nfr.pdf")
dba.plotMA(liverPhe21_peaks, contrast=1, method=DBA_EDGER)
dev.off()

# Volcano plots
dba.plotVolcano(liverPhe21_peaks, contrast=1, method=DBA_DESEQ2)
dba.plotVolcano(liverPhe21_peaks, contrast=1, method=DBA_EDGER)

# Boxplots
sum(liverPhe21_peaks.DB$Fold<0)
sum(liverPhe21_peaks.DB$Fold>0)

liverPhe21_pvals <- dba.plotBox(liverPhe21_peaks)
liverPhe21_pvals

# Heatmaps
liverPhe21_corvals <- dba.plotHeatmap(liverPhe21_peaks)
liverPhe21_hmap <- colorRampPalette(c("red","black","green"))(n=13)
liverPhe21_readscores <- dba.plotHeatmap(liverPhe21_peaks, contrast=1, correlations=FALSE, scale="row", colScheme=liverPhe21_hmap)



### LOESS NORMALISED TEST ###
library(csaw)
loess.liverPhe21_peaks <- liverPhe21_peaks
loess.liverPhe21_peaks$config$AnalysisMethod <- DBA_DESEQ2
loess.liverPhe21_peaks <- dba.normalize(loess.liverPhe21_peaks, offsets=TRUE)
loess.liverPhe21_peaks <- dba.contrast(loess.liverPhe21_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
loess.liverPhe21_peaks <- dba.analyze(loess.liverPhe21_peaks, bGreylist=FALSE)

pdf("loess.liverPhe21_peaks.EDGER.MAplot.nfr.pdf")
dba.plotMA(loess.liverPhe21_peaks)
dev.off()
loess.liverPhe21_peaks.DB.EDGER <- dba.report(loess.liverPhe21_peaks)
sum(loess.liverPhe21_peaks.DB.EDGER$Fold > 0)
sum(loess.liverPhe21_peaks.DB.EDGER$Fold < 0)

pdf("loess.liverPhe21_peaks.DESEQ2.MAplot.nfr.pdf")
dba.plotMA(loess.liverPhe21_peaks)
dev.off()
loess.liverPhe21_peaks.DB.DESEQ2 <- dba.report(loess.liverPhe21_peaks)
sum(loess.liverPhe21_peaks.DB.DESEQ2$Fold > 0)
sum(loess.liverPhe21_peaks.DB.DESEQ2$Fold < 0)

write.table(loess.liverPhe21_peaks.DB.EDGER, "liverPhe21_v_WT.high_conf_peaks.loess_EDGER.DB_results.nfr.txt", quote=FALSE, sep="\t")
write.table(loess.liverPhe21_peaks.DB.DESEQ2, "liverPhe21_v_WT.high_conf_peaks.loess_DESEQ2.DB_results.nfr.txt", quote=FALSE, sep="\t")


writeLines(capture.output(liverPhe21_peaks), "liverPhe21_peaks.nfr.txt")

