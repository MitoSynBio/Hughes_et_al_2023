wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
datadir=$wd/data/AGRF_CAGRF21118859_HVGVVDRXY
trimdir=$wd/pipeline/1_trimming/AGRF_CAGRF21118859_HVGVVDRXY
alndir=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY
index=/media/Data/ref_seqs/Mus_musculus/mm10_GRCm38/GENCODE/bowtie2_GRCm38_GENCODE/bowtie2_GRCm38_GENCODE
picard=/media/Data/tools/picard-tools/picard-2.23.8/picard.jar

## Align (suppressing discordant alignments, unpaired alignments and insert sizes > 1000 are considered discordant; these would be filtered out anyway)
for tissue in liver
do
	for sample in $datadir/*_${tissue}*_L00*_R1.fastq.gz
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}_ATACseq_*_L00*_R1.fastq.gz})"`
		r1=$trimdir/${tissue}/${samplename}/${samplename}_${tissue}_*_R1_val_1.fq.gz
		r2=$trimdir/${tissue}/${samplename}/${samplename}_${tissue}_*_R2_val_2.fq.gz
		output=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.dup.bam
		log=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.log

		mkdir $alndir/${tissue}/${samplename}_${tissue}
		bowtie2 --very-sensitive --dovetail --no-discordant --no-mixed -p 16 -k 20 -X 1000 --fr -x $index -1 $r1 -2 $r2 2> $log | samtools sort -@ 12 -o $output
		samtools index -@ 20 $output

		tmp=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.tmp.bam
		metrics=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.dup.metrics

		java -jar $picard MarkDuplicates I=$output O=$tmp M=$metrics REMOVE_DUPLICATES=FALSE
		rm $output $output.bai
		mv $tmp $output
		samtools index -@ 20 $output
	done
done

## Filtering (blacklist is chrM and list from ENCODE)
for tissue in liver
do
	for sample in $datadir/*_${tissue}*_L00*_R1.fastq.gz
	do
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}_ATACseq_*_L00*_R1.fastq.gz})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.dup.bam
		blackregions=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.blacklisted_regions.bam

		tmp1=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.tmp1.bam
		tmp2=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.tmp2.bam
		output=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.bam
		blacklist=$alndir/blacklist.bed
		metrics=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.metrics

		samtools view -@ 12 -q 30 -f 2 -F 256 -F 4 -o $tmp1 $input
		samtools index $tmp1
		samtools view -@ 12 -L $blacklist -U $tmp2 -o $blackregions $tmp1
		rm $tmp1 $tmp1.bai
		samtools index $tmp2
		samtools index $blackregions

		java -jar $picard MarkDuplicates I=$tmp2 O=$output M=$metrics REMOVE_DUPLICATES=TRUE
		samtools index $output
		rm $tmp2 $tmp2.bai
	done
done

## ATACGraph ##

wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
ag=/media/Data/tools/ATACgraph/ATACgraph/script

## calculate periodicity - round 3
trimdir=$wd/pipeline/1_trimming/AGRF_CAGRF21118859_HVGVVDRXY
alndir=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY
graphdir=$wd/pipeline/ATACgraph/round_3

datadir=$wd/data/AGRF_CAGRF21118859_HVGVVDRXY
for sample in $datadir/*_liver*_L00*_R1.fastq.gz
do
	samplename=`echo "$(b=${sample##*/}; echo ${b%_liver_ATACseq_*_L00*_R1.fastq.gz})"`
	alndir=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY
	input=$alndir/liver/${samplename}_liver/${samplename}_liver.bowtie2.filt.bam
	outdir=$graphdir/liver/${samplename}
	outdist=$outdir/${samplename}_${tissue}.fragdist.pdf
	outfft=$outdir/${samplename}_${tissue}.fragfft.pdf

	mkdir -p $outdir
	$ag/ATACgraph 01_calFragDist $input $outdist $outfft
done


## Filter for nucleosome-free regions (< 175 nt fragment size - determined by ATACgraph)
for tissue in liver
do
	for sample in $datadir/*_${tissue}*_L00*_R1.fastq.gz
	do
		dt=/media/Data/tools/deepTools/deepTools-3.5.0/bin
		samplename=`echo "$(b=${sample##*/}; echo ${b%_${tissue}_ATACseq_*_L00*_R1.fastq.gz})"`
		input=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.filt.bam
		tmp=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.tmp.bam
		output=$alndir/${tissue}/${samplename}_${tissue}/${samplename}_${tissue}.bowtie2.nucleosome_free.175nt.bam

		$dt/alignmentSieve -b $input -o $tmp -p 8 --maxFragmentLength 175
		samtools sort -@ 8 -o $output $tmp
		samtools index $output
		rm $tmp
	done
done



