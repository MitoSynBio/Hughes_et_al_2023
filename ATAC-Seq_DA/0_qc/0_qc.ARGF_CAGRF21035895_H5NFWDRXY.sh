
wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
datadir=$wd/data/AGRF_CAGRF21035895_H5NFWDRXY
qcdir=$wd/pipeline/0_qc

for tissue in brain liver
do
	for sample in $datadir/*_${tissue}*_L00*_R1.fastq.gz
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}_ATACseq_*_L00*_R1.fastq.gz})"`
		r1=$datadir/${samplename}_${tissue}_ATACseq_*_L00*_R1.fastq.gz
		r2=$datadir/${samplename}_${tissue}_ATACseq_*_L00*_R2.fastq.gz
		outdir=$qcdir/round_2/${tissue}/${samplename}_${replicate}

		mkdir $outdir
		fastqc -t 2 -o $outdir $r1 $r2
	done
done

