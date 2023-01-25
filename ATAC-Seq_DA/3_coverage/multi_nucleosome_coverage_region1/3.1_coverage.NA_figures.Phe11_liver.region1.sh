wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
datadir=$wd/data/AGRF_CAGRF21035895_H5NFWDRXY
trimdir=$wd/pipeline/1_trimming
alndir=$wd/pipeline/2_alignment
covdir=$wd/pipeline/3.1_coverage
deeptools=/media/Data/tools/deepTools/deepTools-3.5.0/bin


## Liver Phe11
dir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/3.1_coverage/new_analysis/liverPhe11

WT1=S1_WT
WT2=S3_WT
WT3=22_WT4
WT4=24_WT6
KO1=S4_Phe11
KO2=S6_Phe11
KO3=20_KO5_1_1
KO4=21_KO6_1_1

## NFR
WT1_nfr=$dir/${WT1}_liver.norm_rpm.1_175.extend.atac_shift.bw
WT2_nfr=$dir/${WT2}_liver.norm_rpm.1_175.extend.atac_shift.bw
WT3_nfr=$dir/${WT3}_liver.norm_rpm.1_175.extend.atac_shift.bw
WT4_nfr=$dir/${WT4}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO1_nfr=$dir/${KO1}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO2_nfr=$dir/${KO2}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO3_nfr=$dir/${KO3}_liver.norm_rpm.1_175.extend.atac_shift.bw
KO4_nfr=$dir/${KO4}_liver.norm_rpm.1_175.extend.atac_shift.bw
WT_nfr_out=$dir/WT_liver.norm_rpm.1_175.extend.atac_shift.avg.region1.txt
KO_nfr_out=$dir/KO_liver.norm_rpm.1_175.extend.atac_shift.avg.region1.txt
WTout=$dir/WT_liver.norm_rpm.1_175.extend.atac_shift.avg.region1.npz
KOout=$dir/KO_liver.norm_rpm.1_175.extend.atac_shift.avg.region1.npz

multiBigwigSummary bins -b $WT1_nfr $WT2_nfr $WT3_nfr $WT4_nfr -o $WTout --outRawCounts $WT_nfr_out -p 12 -bs 1 --region chr5:125269989:125526293
multiBigwigSummary bins -b $KO1_nfr $KO2_nfr $KO3_nfr $KO4_nfr -o $KOout --outRawCounts $KO_nfr_out -p 12 -bs 1 --region chr5:125269989:125526293

WTnfr_bg=$dir/WT_liver.norm_rpm.1_175.extend.atac_shift.avg.region1.bedgraph
KOnfr_bg=$dir/KO_liver.norm_rpm.1_175.extend.atac_shift.avg.region1.bedgraph

tail -n+2 $WT_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $WTnfr_bg
tail -n+2 $KO_nfr_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOnfr_bg

nfr_lfc=$dir/Phe11_liver.norm_rpm.1_175.extend.atac_shift.lfc.region1.bedgraph
paste $WTnfr_bg $KOnfr_bg | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3, log(($8+1)/($4+1))/log(2)}' > $nfr_lfc

sed -i '1s/^/track type=bedGraph name="WT_average.nfr.liver"\n/' $WTnfr_bg
sed -i '1s/^/track type=bedGraph name="KO_average.nfr.liver"\n/' $KOnfr_bg
sed -i '1s/^/track type=bedGraph name="log2FoldChange.nfr.liver"\n/' $nfr_lfc



## mono
WT1_mono=$dir/${WT1}_liver.norm_rpm.176_350.extend.atac_shift.bw
WT2_mono=$dir/${WT2}_liver.norm_rpm.176_350.extend.atac_shift.bw
WT3_mono=$dir/${WT3}_liver.norm_rpm.176_350.extend.atac_shift.bw
WT4_mono=$dir/${WT4}_liver.norm_rpm.176_350.extend.atac_shift.bw
KO1_mono=$dir/${KO1}_liver.norm_rpm.176_350.extend.atac_shift.bw
KO2_mono=$dir/${KO2}_liver.norm_rpm.176_350.extend.atac_shift.bw
KO3_mono=$dir/${KO3}_liver.norm_rpm.176_350.extend.atac_shift.bw
KO4_mono=$dir/${KO4}_liver.norm_rpm.176_350.extend.atac_shift.bw
WT_mono_out=$dir/WT_liver.norm_rpm.176_350.extend.atac_shift.avg.region1.txt
KO_mono_out=$dir/KO_liver.norm_rpm.176_350.extend.atac_shift.avg.region1.txt
WTout=$dir/WT_liver.norm_rpm.176_350.extend.atac_shift.avg.region1.npz
KOout=$dir/KO_liver.norm_rpm.176_350.extend.atac_shift.avg.region1.npz

multiBigwigSummary bins -b $WT1_mono $WT2_mono $WT3_mono $WT4_mono -o $WTout --outRawCounts $WT_mono_out -p 12 -bs 1 --region chr5:125269989:125526293
multiBigwigSummary bins -b $KO1_mono $KO2_mono $KO3_mono $KO4_mono -o $KOout --outRawCounts $KO_mono_out -p 12 -bs 1 --region chr5:125269989:125526293

WTmono_bg=$dir/WT_liver.norm_rpm.176_350.extend.atac_shift.avg.region1.bedgraph
KOmono_bg=$dir/KO_liver.norm_rpm.176_350.extend.atac_shift.avg.region1.bedgraph

tail -n+2 $WT_mono_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $WTmono_bg
tail -n+2 $KO_mono_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOmono_bg

mono_lfc=$dir/Phe11_liver.norm_rpm.176_350.extend.atac_shift.lfc.region1.bedgraph
paste $WTmono_bg $KOmono_bg | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3, log(($8+1)/($4+1))/log(2)}' > $mono_lfc

sed -i '1s/^/track type=bedGraph name="WT_average.mono.liver"\n/' $WTmono_bg
sed -i '1s/^/track type=bedGraph name="KO_average.mono.liver"\n/' $KOmono_bg
sed -i '1s/^/track type=bedGraph name="log2FoldChange.mono.liver"\n/' $mono_lfc


## di
WT1_di=$dir/${WT1}_liver.norm_rpm.351_525.extend.atac_shift.bw
WT2_di=$dir/${WT2}_liver.norm_rpm.351_525.extend.atac_shift.bw
WT3_di=$dir/${WT3}_liver.norm_rpm.351_525.extend.atac_shift.bw
WT4_di=$dir/${WT4}_liver.norm_rpm.351_525.extend.atac_shift.bw
KO1_di=$dir/${KO1}_liver.norm_rpm.351_525.extend.atac_shift.bw
KO2_di=$dir/${KO2}_liver.norm_rpm.351_525.extend.atac_shift.bw
KO3_di=$dir/${KO3}_liver.norm_rpm.351_525.extend.atac_shift.bw
KO4_di=$dir/${KO4}_liver.norm_rpm.351_525.extend.atac_shift.bw
WT_di_out=$dir/WT_liver.norm_rpm.351_525.extend.atac_shift.avg.region1.txt
KO_di_out=$dir/KO_liver.norm_rpm.351_525.extend.atac_shift.avg.region1.txt
WTout=$dir/WT_liver.norm_rpm.351_525.extend.atac_shift.avg.region1.npz
KOout=$dir/KO_liver.norm_rpm.351_525.extend.atac_shift.avg.region1.npz

multiBigwigSummary bins -b $WT1_di $WT2_di $WT3_di $WT4_di -o $WTout --outRawCounts $WT_di_out -p 12 -bs 1 --region chr5:125269989:125526293
multiBigwigSummary bins -b $KO1_di $KO2_di $KO3_di $KO4_di -o $KOout --outRawCounts $KO_di_out -p 12 -bs 1 --region chr5:125269989:125526293

WTdi_bg=$dir/WT_liver.norm_rpm.351_525.extend.atac_shift.avg.region1.bedgraph
KOdi_bg=$dir/KO_liver.norm_rpm.351_525.extend.atac_shift.avg.region1.bedgraph

tail -n+2 $WT_di_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $WTdi_bg
tail -n+2 $KO_di_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOdi_bg

di_lfc=$dir/Phe11_liver.norm_rpm.351_525.extend.atac_shift.lfc.region1.bedgraph
paste $WTdi_bg $KOdi_bg | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3, log(($8+1)/($4+1))/log(2)}' > $di_lfc

sed -i '1s/^/track type=bedGraph name="WT_average.di.liver"\n/' $WTdi_bg
sed -i '1s/^/track type=bedGraph name="KO_average.di.liver"\n/' $KOdi_bg
sed -i '1s/^/track type=bedGraph name="log2FoldChange.di.liver"\n/' $di_lfc



## tri
WT1_tri=$dir/${WT1}_liver.norm_rpm.526_700.extend.atac_shift.bw
WT2_tri=$dir/${WT2}_liver.norm_rpm.526_700.extend.atac_shift.bw
WT3_tri=$dir/${WT3}_liver.norm_rpm.526_700.extend.atac_shift.bw
WT4_tri=$dir/${WT4}_liver.norm_rpm.526_700.extend.atac_shift.bw
KO1_tri=$dir/${KO1}_liver.norm_rpm.526_700.extend.atac_shift.bw
KO2_tri=$dir/${KO2}_liver.norm_rpm.526_700.extend.atac_shift.bw
KO3_tri=$dir/${KO3}_liver.norm_rpm.526_700.extend.atac_shift.bw
KO4_tri=$dir/${KO4}_liver.norm_rpm.526_700.extend.atac_shift.bw
WT_tri_out=$dir/WT_liver.norm_rpm.526_700.extend.atac_shift.avg.region1.txt
KO_tri_out=$dir/KO_liver.norm_rpm.526_700.extend.atac_shift.avg.region1.txt
WTout=$dir/WT_liver.norm_rpm.526_700.extend.atac_shift.avg.region1.npz
KOout=$dir/KO_liver.norm_rpm.526_700.extend.atac_shift.avg.region1.npz

multiBigwigSummary bins -b $WT1_tri $WT2_tri $WT3_tri $WT4_tri -o $WTout --outRawCounts $WT_tri_out -p 12 -bs 1 --region chr5:125269989:125526293
multiBigwigSummary bins -b $KO1_tri $KO2_tri $KO3_tri $KO4_tri -o $KOout --outRawCounts $KO_tri_out -p 12 -bs 1 --region chr5:125269989:125526293

WTtri_bg=$dir/WT_liver.norm_rpm.526_700.extend.atac_shift.avg.region1.bedgraph
KOtri_bg=$dir/KO_liver.norm_rpm.526_700.extend.atac_shift.avg.region1.bedgraph

tail -n+2 $WT_tri_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $WTtri_bg
tail -n+2 $KO_tri_out | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3,($4+$5+$6+$7)/4}' > $KOtri_bg

tri_lfc=$dir/Phe11_liver.norm_rpm.526_700.extend.atac_shift.lfc.region1.bedgraph
paste $WTtri_bg $KOtri_bg | awk 'BEGIN {FS="\t";OFS="\t"}; {print $1,$2,$3, log(($8+1)/($4+1))/log(2)}' > $tri_lfc

sed -i '1s/^/track type=bedGraph name="WT_average.tri.liver"\n/' $WTtri_bg
sed -i '1s/^/track type=bedGraph name="KO_average.tri.liver"\n/' $KOtri_bg
sed -i '1s/^/track type=bedGraph name="log2FoldChange.tri.liver"\n/' $tri_lfc


