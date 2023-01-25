
########### round 2 
wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
datadir=$wd/data/AGRF_CAGRF21035895_H5NFWDRXY
trimdir=$wd/pipeline/1_trimming/round_2
alndir=$wd/pipeline/2_alignment/round_2
covdir=$wd/pipeline/3.1_coverage/round_2
deeptools=/media/Data/tools/deepTools/deepTools-3.5.0/bin

for tissue in liver
do
	for sample in $alndir/${tissue}/*/*_${tissue}*.bowtie2.filt.bam
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}.bowtie2.filt.bam})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.bam
		output=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam

		$deeptools/alignmentSieve -p 16 -b $input -o /dev/stdout --ATACshift | samtools sort -@ 16 -O BAM -o $output
		samtools index $output
	done
done


# normalised coverage - use filt
for tissue in brain liver
do
	for sample in $alndir/${tissue}/*/*_${tissue}*.bowtie2.filt.bam
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}.bowtie2.filt.bam})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.atac_shift.bw

		mkdir -p $outdir
		$deeptools/bamCoverage --normalizeUsing CPM -p 10 -bs 1 -of bigwig -b $input -o $output
	done
done


#############################################################
# nucleosome-free (175nt) - extend reads (so that coverage shows nucleosome-free and -bound regions)
for tissue in liver
do
	for sample in $alndir/${tissue}/*/*_${tissue}*.bowtie2.filt.bam
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}.bowtie2.filt.bam})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.1_175.extend.atac_shift.bw

		$deeptools/bamCoverage --extendReads --minFragmentLength 1 --maxFragmentLength 175 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done

# nucleosome-bound (176nt +) - extend reads (so that coverage shows nucleosome-free and -bound regions)
for tissue in liver
do
	for sample in $alndir/${tissue}/*/*_${tissue}*.bowtie2.filt.bam
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}.bowtie2.filt.bam})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.176_above.extend.atac_shift.bw

		$deeptools/bamCoverage --extendReads --minFragmentLength 176 --maxFragmentLength 1000 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done

## 176-350 (mononucleosome-bound) - extend
for tissue in liver
do
	for sample in $alndir/${tissue}/*/*_${tissue}*.bowtie2.filt.bam
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}.bowtie2.filt.bam})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.176_350.extend.atac_shift.bw

		$deeptools/bamCoverage --extendReads --minFragmentLength 176 --maxFragmentLength 350 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done

## 351-525 (dinucleosome-bound) - extend
for tissue in liver
do
	for sample in $alndir/${tissue}/*/*_${tissue}*.bowtie2.filt.bam
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}.bowtie2.filt.bam})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.351_525.extend.atac_shift.bw

		bamCoverage --extendReads --minFragmentLength 351 --maxFragmentLength 525 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done

## 526-700 (trinucleosome-bound)
for tissue in brain
do
	for sample in $alndir/${tissue}/*/*_${tissue}*.bowtie2.filt.bam
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}.bowtie2.filt.bam})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.526_700.extend.atac_shift.bw

		$deeptools/bamCoverage --extendReads --minFragmentLength 526 --maxFragmentLength 700 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done


## 701-875 (tetranucleosome-bound)
for tissue in liver
do
	for sample in $alndir/${tissue}/*/*_${tissue}*.bowtie2.filt.bam
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}.bowtie2.filt.bam})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.701_875.extend.atac_shift.bw

		$deeptools/bamCoverage --extendReads --minFragmentLength 701 --maxFragmentLength 875 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done



