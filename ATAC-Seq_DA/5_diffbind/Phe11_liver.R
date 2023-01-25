
setwd("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis")

library(DiffBind)
library(rtracklayer)

samples <- read.table("/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/5_diffBind/R_4.0.3_hc/canon/new_analysis/Phe11_WT.liver.samples.txt",header=TRUE,sep="\t")

#### liverPhe11 SAMPLE ANALYSIS - will aid in determining quality of results - PCA all vs BD - account for tissue ####
## 1. Read in liverPhe11_peaksets
liverPhe11_peaks <- dba(sampleSheet=samples)
liverPhe11_peaks

plot(liverPhe11_peaks)

## 1a. Blacklists (and greylist?)
liverPhe11_peakdata <- dba.show(liverPhe11_peaks)$Intervals
liverPhe11_peaks <- dba.blacklist(liverPhe11_peaks, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE)
liverPhe11_peaks

liverPhe11_peakdata.BL <- dba.show(liverPhe11_peaks)$Intervals
liverPhe11_peakdata - liverPhe11_peakdata.BL

## 1b. Consensus peaks
liverPhe11_consensus <- dba.peakset(liverPhe11_peaks, consensus=DBA_TISSUE, minOverlap=1)	# for replicate consensus, but has already been done by genrich
liverPhe11_consensus <- dba(liverPhe11_consensus, mask=liverPhe11_consensus$masks$"Consensus", minOverlap=1)
liverPhe11_consensus_peaks <- dba.peakset(liverPhe11_consensus, bRetrieve=TRUE)
liverPhe11_consensus_peaks

write.table(liverPhe11_consensus_peaks, "liverPhe11_v_WT_consensus_peaks.high_conf_peaks.nfr.txt", sep="\t")

## 2. Counting reads
liverPhe11_peaks <- dba.count(liverPhe11_peaks, peaks=liverPhe11_consensus_peaks, bParallel=FALSE, summits=FALSE, mapQCth=30, fragmentSize=0)
liverPhe11_peaks

liverPhe11_info <- dba.show(liverPhe11_peaks)
liverPhe11_libsizes <- cbind(LibReads=liverPhe11_info$Reads, FRiP=liverPhe11_info$FRiP, PeakReads=round(liverPhe11_info$Reads * liverPhe11_info$FRiP))
rownames(liverPhe11_libsizes) <- liverPhe11_info$ID
liverPhe11_libsizes

pdf("liverPhe11/liverPhe11_peaks.high_conf_peaks.all_sites.nfr.pdf")
plot(liverPhe11_peaks)
dev.off()

## 3. Normalisation (try all three)
liverPhe11_peaks <- dba.normalize(liverPhe11_peaks, method=DBA_ALL_METHODS, normalize=DBA_NORM_NATIVE) ## by default library is set to RiP
liverPhe11_norm.DESeq2 <- dba.normalize(liverPhe11_peaks, bRetrieve=TRUE)
liverPhe11_norm.DESeq2

liverPhe11_norm.edgeR <- dba.normalize(liverPhe11_peaks, bRetrieve=TRUE, method=DBA_EDGER)
liverPhe11_norm.edgeR

liverPhe11_normlibs.DESeq2 <- cbind(FullLibSize=liverPhe11_norm.DESeq2$lib.sizes, NormFacs=liverPhe11_norm.DESeq2$norm.factors, NormLibSize=round(liverPhe11_norm.DESeq2$lib.sizes/liverPhe11_norm.DESeq2$norm.factors))
rownames(liverPhe11_normlibs.DESeq2) <- liverPhe11_info$ID
liverPhe11_normlibs.DESeq2

liverPhe11_normlibs.edgeR <- cbind(FullLibSize=liverPhe11_norm.edgeR$lib.sizes, NormFacs=liverPhe11_norm.edgeR$norm.factors, NormLibSize=round(liverPhe11_norm.edgeR$lib.sizes/liverPhe11_norm.edgeR$norm.factors))
rownames(liverPhe11_normlibs.edgeR) <- liverPhe11_info$ID
liverPhe11_normlibs.edgeR

write.table(liverPhe11_normlibs.DESeq2, "liverPhe11_normlibs_DESeq2.RiP.nfr.txt", sep="\t", quote=FALSE)
write.table(liverPhe11_normlibs.edgeR, "liverPhe11_normlibs_edgeR.RiP.nfr.txt", sep="\t", quote=FALSE)

## 4. Establishing a model design
liverPhe11_peaks <- dba.contrast(liverPhe11_peaks, design="~Condition", reorderMeta=list(Condition="WT"))
liverPhe11_peaks

## 5. Performing the analysis
liverPhe11_peaks <- dba.analyze(liverPhe11_peaks, method=DBA_ALL_METHODS, bGreylist=FALSE)
dba.show(liverPhe11_peaks, bContrasts=TRUE)
liverPhe11_results_DESeq2 <- dba.analyze(liverPhe11_peaks, bRetrieveAnalysis=DBA_DESEQ2)
liverPhe11_results_edgeR <- dba.analyze(liverPhe11_peaks, bRetrieveAnalysis=DBA_EDGER)


pdf("liverPhe11/liverPhe11_peaks.high_conf_peaks.DESeq2.DB_peaks.nfr.pdf")
plot(liverPhe11_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("liverPhe11/liverPhe11_peaks.high_conf_peaks.edgeR.DB_peaks.nfr.pdf")
plot(liverPhe11_peaks, contrast=1, method=DBA_EDGER)
dev.off()

## 6. Retrieving the DB sites
liverPhe11_peaks.DB.DESeq2 <- dba.report(liverPhe11_peaks, contrast=1, method=DBA_DESEQ2)
sum(liverPhe11_peaks.DB.DESeq2$Fold>0)
sum(liverPhe11_peaks.DB.DESeq2$Fold<0)

liverPhe11_peaks.DB.edgeR <- dba.report(liverPhe11_peaks, contrast=1, method=DBA_EDGER)
sum(liverPhe11_peaks.DB.edgeR$Fold>0)
sum(liverPhe11_peaks.DB.edgeR$Fold<0)

write.table(liverPhe11_peaks.DB.DESeq2, "liverPhe11_v_WT.high_conf_peaks.DESEQ2_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")
write.table(liverPhe11_peaks.DB.edgeR, "liverPhe11_v_WT.high_conf_peaks.EDGER_RiP.DB_results.nfr.txt", quote=FALSE, sep="\t")


liverPhe11_peaks.ALL.DESeq2 <- dba.report(liverPhe11_peaks, method=DBA_DESEQ2, th=1)
liverPhe11_peaks.ALL.edgeR <- dba.report(liverPhe11_peaks, method=DBA_EDGER, th=1)
write.table(liverPhe11_peaks.ALL.DESeq2, "liverPhe11_v_WT.high_conf_peaks.DESEQ2_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")
write.table(liverPhe11_peaks.ALL.edgeR, "liverPhe11_v_WT.high_conf_peaks.EDGER_RiP.ALL_results.nfr.txt", quote=FALSE, sep="\t")

save.image("liverPhe11/liverPhe11_v_WT.high_conf_peaks.RData")

writeLines(capture.output(sessionInfo()), "sessionInfo.nfr.txt")

##############################
##			PLOTS			##
##############################

# Venn diagrams
dba.plotVenn(liverPhe11_peaks, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)

# PCA Plots
pdf("liverPhe11/liverPhe11_v_WT.high_conf_peaks.ALLSITES.PCA.nfr.pdf")
dba.plotPCA(liverPhe11_peaks, DBA_TISSUE, label=DBA_CONDITION)
dev.off()
pdf("liverPhe11/liverPhe11_v_WT.high_conf_peaks.DESEQ2.PCA.nfr.pdf")
dba.plotPCA(liverPhe11_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("liverPhe11/liverPhe11_v_WT.high_conf_peaks.EDGER.PCA.nfr.pdf")
dba.plotPCA(liverPhe11_peaks, DBA_TISSUE, label=DBA_CONDITION, contrast=1, method=DBA_EDGER)
dev.off()


#dba.plotPCA(liverPhe11_peaks, contrast=1, label=DBA_TISSUE)
#dba.plotPCA(liverPhe11_peaks, attributes=c(DBA_TISSUE,DBA_CONDITION), label=DBA_REPLICATE)
#dba.plotPCA(liverPhe11_peaks, contrast=1)

# MA plots
pdf("liverPhe11/liverPhe11_v_WT.high_conf_peaks.DESEQ2.MAplot.nfr.pdf")
dba.plotMA(liverPhe11_peaks, contrast=1, method=DBA_DESEQ2)
dev.off()
pdf("liverPhe11/liverPhe11_v_WT.high_conf_peaks.EDGER.MAplot.nfr.pdf")
dba.plotMA(liverPhe11_peaks, contrast=1, method=DBA_EDGER)
dev.off()

# Volcano plots
dba.plotVolcano(liverPhe11_peaks, contrast=1)

# Boxplots
sum(liverPhe11_peaks.DB$Fold<0)
sum(liverPhe11_peaks.DB$Fold>0)

liverPhe11_pvals <- dba.plotBox(liverPhe11_peaks)
liverPhe11_pvals

# Heatmaps
liverPhe11_corvals <- dba.plotHeatmap(liverPhe11_peaks)
liverPhe11_hmap <- colorRampPalette(c("red","black","green"))(n=13)
liverPhe11_readscores <- dba.plotHeatmap(liverPhe11_peaks, contrast=1, correlations=FALSE, scale="row", colScheme=liverPhe11_hmap)


writeLines(capture.output(liverPhe11_peaks), "liverPhe11_peaks.nfr.txt")

