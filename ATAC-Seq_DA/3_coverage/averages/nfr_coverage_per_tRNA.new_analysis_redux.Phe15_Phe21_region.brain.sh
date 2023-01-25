wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
datadir=$wd/data/AGRF_CAGRF21035895_H5NFWDRXY
trimdir=$wd/pipeline/1_trimming
alndir=$wd/pipeline/2_alignment
covdir=$wd/pipeline/3.1_coverage
deeptools=/media/Data/tools/deepTools/deepTools-3.5.0/bin


#### BRAIN ####

round2_brain=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/3.1_coverage/round_2/filt/brain
round1_brain=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/3.1_coverage/filt/brain
outdir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/3.1_coverage/new_analysis/averages/brain


################################################################################################################################
# tRNA-Phe-1-5 & Phe-2-1 gene region #
region="chr19:11987618:12019348"
################################################################################################################################

########
## WT ##
########

WT1=S1_WT
WT2=S2_WT
WT3=S3_WT
WT4=16_WT4
WT5=17_WT5
WT6=18_WT6

## NFR
WT1_nfr=$round1_brain/${WT1}_brain.norm_rpm.1_175.extend.atac_shift.bw
WT2_nfr=$round1_brain/${WT2}_brain.norm_rpm.1_175.extend.atac_shift.bw
WT3_nfr=$round1_brain/${WT3}_brain.norm_rpm.1_175.extend.atac_shift.bw
WT4_nfr=$round2_brain/${WT4}_brain.norm_rpm.1_175.extend.atac_shift.bw
WT5_nfr=$round2_brain/${WT5}_brain.norm_rpm.1_175.extend.atac_shift.bw
WT6_nfr=$round2_brain/${WT6}_brain.norm_rpm.1_175.extend.atac_shift.bw

WT_nfr_out=$outdir/WT_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.txt
WTout=$outdir/WT_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.npz

multiBigwigSummary bins -b $WT1_nfr $WT2_nfr $WT3_nfr $WT4_nfr $WT5_nfr $WT6_nfr -o $WTout --outRawCounts $WT_nfr_out -p 12 -bs 1 --region $region

WTnfr_bg=$outdir/WT_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.bedgraph

tail -n+2 $WT_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $WTnfr_bg
sed -i '1s/^/track type=bedGraph name="WT_average.nfr.brain"\n/' $WTnfr_bg


########
## Phe11 ##
########

KO1=S5_Phe11
KO2=S6_Phe11
KO3=13_KO4_1_1
KO4=14_KO5_1_1
KO5=15_KO6_1_1

## NFR
KO1_nfr=$round1_brain/${KO1}_brain.norm_rpm.1_175.extend.atac_shift.bw
KO2_nfr=$round1_brain/${KO2}_brain.norm_rpm.1_175.extend.atac_shift.bw
KO3_nfr=$round2_brain/${KO3}_brain.norm_rpm.1_175.extend.atac_shift.bw
KO4_nfr=$round2_brain/${KO4}_brain.norm_rpm.1_175.extend.atac_shift.bw
KO5_nfr=$round2_brain/${KO5}_brain.norm_rpm.1_175.extend.atac_shift.bw

KO_nfr_out=$outdir/Phe11_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.txt
KOout=$outdir/Phe11_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.npz

multiBigwigSummary bins -b $KO1_nfr $KO2_nfr $KO3_nfr $KO4_nfr $KO5_nfr -o $KOout --outRawCounts $KO_nfr_out -p 12 -bs 1 --region $region

KOnfr_bg=$outdir/Phe11_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.bedgraph

tail -n+2 $KO_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOnfr_bg
sed -i '1s/^/track type=bedGraph name="Phe11_average.nfr.brain"\n/' $KOnfr_bg



########
## Phe21 ##
########

KO1=S7_Phe21
KO2=S8_Phe21
KO3=S9_Phe21

## NFR
KO1_nfr=$round1_brain/${KO1}_brain.norm_rpm.1_175.extend.atac_shift.bw
KO2_nfr=$round1_brain/${KO2}_brain.norm_rpm.1_175.extend.atac_shift.bw
KO3_nfr=$round1_brain/${KO3}_brain.norm_rpm.1_175.extend.atac_shift.bw

KO_nfr_out=$outdir/Phe21_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.txt
KOout=$outdir/Phe21_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.npz

multiBigwigSummary bins -b $KO1_nfr $KO2_nfr $KO3_nfr -o $KOout --outRawCounts $KO_nfr_out -p 12 -bs 1 --region $region

KOnfr_bg=$outdir/Phe21_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.bedgraph

tail -n+2 $KO_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOnfr_bg
sed -i '1s/^/track type=bedGraph name="Phe21_average.nfr.brain"\n/' $KOnfr_bg




########
## Phe31 ##
########

KO1=S10_Phe31
KO2=S11_Phe31
KO3=S12_Phe31

## NFR
KO1_nfr=$round1_brain/${KO1}_brain.norm_rpm.1_175.extend.atac_shift.bw
KO2_nfr=$round1_brain/${KO2}_brain.norm_rpm.1_175.extend.atac_shift.bw
KO3_nfr=$round1_brain/${KO3}_brain.norm_rpm.1_175.extend.atac_shift.bw

KO_nfr_out=$outdir/Phe31_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.txt
KOout=$outdir/Phe31_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.npz

multiBigwigSummary bins -b $KO1_nfr $KO2_nfr $KO3_nfr -o $KOout --outRawCounts $KO_nfr_out -p 12 -bs 1 --region $region

KOnfr_bg=$outdir/Phe31_brain.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe15_Phe21.bedgraph

tail -n+2 $KO_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOnfr_bg
sed -i '1s/^/track type=bedGraph name="Phe31_average.nfr.brain"\n/' $KOnfr_bg




