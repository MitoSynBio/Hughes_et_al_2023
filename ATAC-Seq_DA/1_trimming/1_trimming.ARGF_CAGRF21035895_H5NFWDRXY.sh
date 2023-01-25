
wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
datadir=$wd/data/AGRF_CAGRF21035895_H5NFWDRXY
trimdir=$wd/pipeline/1_trimming
tg=/media/Data/tools/TrimGalore/TrimGalore-0.6.5/trim_galore

for tissue in brain liver
do
	for sample in $datadir/*_${tissue}*_L00*_R1.fastq.gz
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}_ATACseq_*_L00*_R1.fastq.gz})"`
		r1=$datadir/${samplename}_${tissue}_ATACseq_*_L00*_R1.fastq.gz
		r2=$datadir/${samplename}_${tissue}_ATACseq_*_L00*_R2.fastq.gz
		outdir=$trimdir/round_2/${tissue}/${samplename}

		mkdir $outdir
		$tg --paired --fastqc --output_dir $outdir $r1 $r2
	done
done

