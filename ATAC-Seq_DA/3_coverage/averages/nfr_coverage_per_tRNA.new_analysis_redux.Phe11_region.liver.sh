wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
datadir=$wd/data/AGRF_CAGRF21035895_H5NFWDRXY
trimdir=$wd/pipeline/1_trimming
alndir=$wd/pipeline/2_alignment
covdir=$wd/pipeline/3.1_coverage
deeptools=/media/Data/tools/deepTools/deepTools-3.5.0/bin


#### liver ####

round2_brain=/media/IRDS/Server/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/3.1_coverage/round_2/filt/liver
round1_brain=/media/IRDS/Server/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/3.1_coverage/filt/liver
outdir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/3.1_coverage/new_analysis.redux/averages/liver


################################################################################################################################
# tRNA-Phe-1-1 gene region #
region="chr5:125389790:125417802"
################################################################################################################################


### Add new liver WT samples from new sequencing
########
## WT ##
########

WT1=S1_WT
WT2=S3_WT
WT3=22_WT4
WT4=24_WT6

## NFR
WT1_nfr=$round1_liver/${WT1}_liver.norm_rpm.1_175.extend.atac_shift.bw
WT2_nfr=$round1_liver/${WT2}_liver.norm_rpm.1_175.extend.atac_shift.bw
WT3_nfr=$round2_liver/${WT3}_liver.norm_rpm.1_175.extend.atac_shift.bw
WT4_nfr=$round2_liver/${WT4}_liver.norm_rpm.1_175.extend.atac_shift.bw

WT_nfr_out=$outdir/WT_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.txt
WTout=$outdir/WT_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.npz

multiBigwigSummary bins -b $WT1_nfr $WT2_nfr $WT3_nfr $WT4_nfr -o $WTout --outRawCounts $WT_nfr_out -p 12 -bs 1 --region $region

WTnfr_bg=$outdir/WT_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.bedgraph

tail -n+2 $WT_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $WTnfr_bg
sed -i '1s/^/track type=bedGraph name="WT_average.nfr.liver"\n/' $WTnfr_bg


########
## Phe11 ##
########

KO1=S4_Phe11
KO3=S6_Phe11
KO5=20_KO5_1_1
KO6=21_KO6_1_1

## NFR
KO1_nfr=$round1_liver/${KO1}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO3_nfr=$round1_liver/${KO3}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO5_nfr=$round2_liver/${KO5}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO6_nfr=$round2_liver/${KO6}_liver.norm_rpm.1_175.extend.atac_shift.bw

KO_nfr_out=$outdir/Phe11_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.txt
KOout=$outdir/Phe11_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.npz

multiBigwigSummary bins -b $KO1_nfr $KO3_nfr $KO5_nfr $KO6_nfr -o $KOout --outRawCounts $KO_nfr_out -p 12 -bs 1 --region $region

KOnfr_bg=$outdir/Phe11_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.bedgraph

tail -n+2 $KO_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOnfr_bg
sed -i '1s/^/track type=bedGraph name="Phe11_average.nfr.liver"\n/' $KOnfr_bg



########
## Phe21 ##
########

KO1=S7_Phe21
KO2=S8_Phe21
KO3=S9_Phe21

## NFR
KO1_nfr=$round1_liver/${KO1}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO2_nfr=$round1_liver/${KO2}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO3_nfr=$round1_liver/${KO3}_liver.norm_rpm.1_175.extend.atac_shift.bw

KO_nfr_out=$outdir/Phe21_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.txt
KOout=$outdir/Phe21_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.npz

multiBigwigSummary bins -b $KO1_nfr $KO2_nfr $KO3_nfr -o $KOout --outRawCounts $KO_nfr_out -p 12 -bs 1 --region $region

KOnfr_bg=$outdir/Phe21_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.bedgraph

tail -n+2 $KO_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOnfr_bg
sed -i '1s/^/track type=bedGraph name="Phe21_average.nfr.liver"\n/' $KOnfr_bg



############## this must be replaced with new data
########
## Phe31 ##
########

round3=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY
KO1=10_31_50_1
KO2=11_31_50_2
KO3=12_31_50_3

## NFR
KO1_nfr=$round3/${KO1}_liver/${KO1}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO2_nfr=$round3/${KO2}_liver/${KO2}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO3_nfr=$round3/${KO3}_liver/${KO3}_liver.norm_rpm.1_175.extend.atac_shift.bw

KO_nfr_out=$outdir/Phe31_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.txt
KOout=$outdir/Phe31_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.npz

multiBigwigSummary bins -b $KO1_nfr $KO2_nfr $KO3_nfr -o $KOout --outRawCounts $KO_nfr_out -p 12 -bs 1 --region $region

KOnfr_bg=$outdir/Phe31_liver.norm_rpm.NFR.extend.atac_shift.avg.tRNAPhe11.bedgraph

tail -n+2 $KO_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOnfr_bg
sed -i '1s/^/track type=bedGraph name="Phe31_average.nfr.liver"\n/' $KOnfr_bg




