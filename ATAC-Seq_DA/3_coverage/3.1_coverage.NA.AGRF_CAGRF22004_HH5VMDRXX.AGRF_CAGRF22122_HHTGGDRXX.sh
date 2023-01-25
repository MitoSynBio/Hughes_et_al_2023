wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
datadir=$wd/data/AGRF_CAGRF21035895_H5NFWDRXY
trimdir=$wd/pipeline/1_trimming
alndir=$wd/pipeline/2_alignment
covdir=$wd/pipeline/3.1_coverage
deeptools=/media/Data/tools/deepTools/deepTools-3.5.0/bin

### round 1
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.bam
		samtools index $input
	done
done

for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.bam
		output=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam

		$deeptools/alignmentSieve -p 16 -b $input -o /dev/stdout --ATACshift | samtools sort -@ 16 -O BAM -o $output
		samtools index $output
	done
done

# normalised coverage - use filt
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.atac_shift.atac_shift.bw

		mkdir -p $outdir
		$deeptools/bamCoverage --normalizeUsing CPM -p 10 -bs 1 -of bigwig -b $input -o $output
	done
done

# nucleosome-free (<= 175nt) - extend reads (so that coverage shows nucleosome-free and -bound regions)
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.1_175.extend.atac_shift.bw

		mkdir -p $outdir
		$deeptools/bamCoverage --extendReads --maxFragmentLength 175 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done

# nucleosome-bound (176nt +) - extend reads (so that coverage shows nucleosome-free and -bound regions)
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.176_above.extend.atac_shift.bw

		mkdir -p $outdir
		$deeptools/bamCoverage --extendReads --minFragmentLength 176 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done

## 176-350 (mononucleosome-bound) - extend
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.176_350.extend.atac_shift.bw

		mkdir -p $outdir
		$deeptools/bamCoverage --extendReads --minFragmentLength 176 --maxFragmentLength 350 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done

## 351-525 (dinucleosome-bound) - extend
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.351_525.extend.atac_shift.bw

		mkdir -p $outdir
		bamCoverage --extendReads --minFragmentLength 351 --maxFragmentLength 525 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done

## 526-700 (trinucleosome-bound)
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.526_700.extend.atac_shift.bw

		mkdir -p $outdir
		$deeptools/bamCoverage --extendReads --minFragmentLength 526 --maxFragmentLength 700 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done


## 701-875 (tetranucleosome-bound)
for tissue in brain liver
do
	for sample in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.atac_shift.bam
		outdir=$covdir/filt/${tissue}

		# library size normalisation (CPM)
		output=$outdir/${samplename}_${tissue}.norm_rpm.701_875.extend.atac_shift.bw

		mkdir -p $outdir
		$deeptools/bamCoverage --extendReads --minFragmentLength 701 --maxFragmentLength 875 --normalizeUsing CPM -p 16 -bs 1 -of bigwig -b $input -o $output
	done
done



