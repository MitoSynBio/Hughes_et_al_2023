WD=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
DATA=$WD/pipeline/1_trimming
INDEX=/media/Data/ref_seqs/Mus_musculus/mm10_GRCm38/GENCODE/bowtie2_GRCm38_GENCODE/bowtie2_GRCm38_GENCODE
PICARD=/media/Data/tools/picard-tools/picard-2.23.8/picard.jar

## Align (suppressing discordant alignments, unpaired alignments and insert sizes > 1000 are considered discordant; these would be filtered out anyway)
for SAMPLE in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
do
	for TISSUE in brain liver
	do
		L001_R1=$DATA/${TISSUE}/${SAMPLE}_L001/${SAMPLE}_${TISSUE}_L001_R1_val_1.fq.gz
		L001_R2=$DATA/${TISSUE}/${SAMPLE}_L001/${SAMPLE}_${TISSUE}_L001_R2_val_2.fq.gz
		L002_R1=$DATA/${TISSUE}/${SAMPLE}_L002/${SAMPLE}_${TISSUE}_L002_R1_val_1.fq.gz
		L002_R2=$DATA/${TISSUE}/${SAMPLE}_L002/${SAMPLE}_${TISSUE}_L002_R2_val_2.fq.gz
		OUTPUT=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.dup.bam
		LOG=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.log

		mkdir $WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}
		bowtie2 --very-sensitive --dovetail --no-discordant --no-mixed -p 16 -k 20 -X 1000 --fr -x $INDEX -1 $L001_R1,$L002_R1 -2 $L001_R2,$L002_R2 2> $LOG | samtools sort -@ 12 -o $OUTPUT
		samtools index -@ 20 $OUTPUT
	done
done

## Mark duplicates
for SAMPLE in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
do
	for TISSUE in brain liver
	do
		input=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.dup.bam
		tmp=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.tmp.bam
		metrics=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.dup.metrics

		java -jar $PICARD MarkDuplicates I=$input O=$tmp M=$metrics REMOVE_DUPLICATES=FALSE
		rm $input $input.bai
		mv $tmp $input
		samtools index -@ 20 $input
	done
done


## Filtering (blacklist is chrM and list from ENCODE)
for TISSUE in brain liver
do
	for SAMPLE in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		INPUT=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.dup.bam
		BLACKREGIONS=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.blacklisted_regions.bam

		TMP1=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.tmp1.bam
		TMP2=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.tmp2.bam
		OUTPUT=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.filt.bam
		BLACKLIST=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/2_alignment/blacklist.bed
		METRICS=$WD/pipeline/2_alignment/${TISSUE}/${SAMPLE}_${TISSUE}/${SAMPLE}_${TISSUE}.bowtie2.filt.metrics

		samtools view -@ 12 -q 30 -f 2 -F 256 -F 4 -o $TMP1 $INPUT
		samtools index $TMP1
		samtools view -@ 12 -L $BLACKLIST -U $TMP2 -o $BLACKREGIONS $TMP1
		rm $TMP1 $TMP1.bai
		samtools index $TMP2
		samtools index $BLACKREGIONS

		java -jar $PICARD MarkDuplicates I=$TMP2 O=$OUTPUT M=$METRICS REMOVE_DUPLICATES=TRUE
		samtools index $OUTPUT
		rm $TMP2 $TMP2.bai
	done
done

## ATACGraph ##
datadir=$wd/data/links
trimdir=$wd/pipeline/1_trimming
alndir=$wd/pipeline/2_alignment
graphdir=$wd/pipeline/ATACgraph/round_1

wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
ag=/media/Data/tools/ATACgraph/ATACgraph/script


## calculate periodicity - round 1
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		alndir=$wd/pipeline/2_alignment
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.bam
		outdir=$graphdir/${tissue}/${samplename}
		outdist=$outdir/${samplename}_${tissue}.fragdist.pdf
		outfft=$outdir/${samplename}_${tissue}.fragfft.pdf

		mkdir -p $outdir
		$ag/ATACgraph 01_calFragDist $input $outdist $outfft
	done
done

## Filter for nucleosome-free regions (< 175 nt fragment size - determined by ATACgraph)
for tissue in brain liver
do
	for samplename in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
	do
		dt=/media/Data/tools/deepTools/deepTools-3.5.0/bin
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.bam
		tmp=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.tmp.bam
		output=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.nucleosome_free.175nt.bam

		$dt/alignmentSieve -b $input -o $tmp -p 8 --maxFragmentLength 175
		samtools sort -@ 8 -o $output $tmp
		samtools index $output
		rm $tmp
	done
done






