wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
genrich=/media/Data/tools/Genrich/Genrich

# Genrich calling peaks for each sample group (genrich is left to do the filtering; manually filtered bams used for coverage profiling)
# blacklisting performed in diffbind downstream
# (-q) qvalue cutoff set at 0.05, (-m) mapping quality required phred 30, (-r) remove PCR duplicates, -(-a) min AUC for a peak =200 (default)
# (-y) keep unpaired alignments (phread quality will ensure good mapping)
## WT for Phe11
for tissue in liver
do
	#WT
	blacklist=$wd/pipeline/2_alignment/ENCFF547MET.bed
	outdir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/4_peakCalling/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}

	WT_R1=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/7_WT_50_1_${tissue}/7_WT_50_1_${tissue}.bowtie2.nucleosome_free.175nt.bam
	WT_R2=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/8_WT_50_2_${tissue}/8_WT_50_2_${tissue}.bowtie2.nucleosome_free.175nt.bam
	WT_R3=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/9_WT_50_3_${tissue}/9_WT_50_3_${tissue}.bowtie2.nucleosome_free.175nt.bam
	R1_tmp=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/7_WT_50_1_${tissue}/7_WT_50_1_${tissue}.bowtie2.tmp.bam
	R2_tmp=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/8_WT_50_2_${tissue}/8_WT_50_2_${tissue}.bowtie2.tmp.bam
	R3_tmp=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/9_WT_50_3_${tissue}/9_WT_50_3_${tissue}.bowtie2.tmp.bam

	WT_peaks=$outdir/WT_${tissue}.forPhe31.genrich.hc.final.nfr.narrowPeak
	WT_log=$outdir/WT_${tissue}.forPhe31.genrich.hc.final.nfr.log
	WT_pu=$outdir/WT_${tissue}.forPhe31.genrich.pileup.hc.final.nfr.bedgraph

	samtools sort -@ 12 -n -O BAM -o $R1_tmp $WT_R1
	samtools sort -@ 12 -n -O BAM -o $R2_tmp $WT_R2
	samtools sort -@ 12 -n -O BAM -o $R3_tmp $WT_R3

	$genrich -t $R1_tmp,$R2_tmp,$R3_tmp -o $WT_peaks -E $blacklist -k $WT_pu -a 200 -m 30 -q 0.05 -j -r -e chrM,chrY -v 2> $WT_log
	gzip $WT_pu

	rm $R1_tmp $R2_tmp $R3_tmp
done

### Phe31
for tissue in liver
do
	#Phe31
	blacklist=$wd/pipeline/2_alignment/ENCFF547MET.bed
	outdir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/4_peakCalling/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}

	Phe31_R1=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/10_31_50_1_${tissue}/10_31_50_1_${tissue}.bowtie2.nucleosome_free.175nt.bam
	Phe31_R2=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/11_31_50_2_${tissue}/11_31_50_2_${tissue}.bowtie2.nucleosome_free.175nt.bam
	Phe31_R3=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/12_31_50_3_${tissue}/12_31_50_3_${tissue}.bowtie2.nucleosome_free.175nt.bam
	R1_tmp=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/10_31_50_1_${tissue}/10_31_50_1_${tissue}.bowtie2.tmp.bam
	R2_tmp=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/11_31_50_2_${tissue}/11_31_50_2_${tissue}.bowtie2.tmp.bam
	R3_tmp=$wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/${tissue}/12_31_50_3_${tissue}/12_31_50_3_${tissue}.bowtie2.tmp.bam

	Phe31_peaks=$outdir/Phe31_${tissue}.genrich.hc.final.nfr.narrowPeak
	Phe31_log=$outdir/Phe31_${tissue}.genrich.hc.final.nfr.log
	Phe31_pu=$outdir/Phe31_${tissue}.genrich.pileup.hc.final.nfr.bedgraph

	samtools sort -@ 12 -n -O BAM -o $R1_tmp $Phe31_R1
	samtools sort -@ 12 -n -O BAM -o $R2_tmp $Phe31_R2
	samtools sort -@ 12 -n -O BAM -o $R3_tmp $Phe31_R3

	$genrich -t $R1_tmp,$R2_tmp,$R3_tmp -o $Phe31_peaks -E $blacklist -k $Phe31_pu -a 200 -m 30 -q 0.05 -j -r -e chrM,chrY -v 2> $Phe31_log
	gzip $Phe31_pu

	rm $R1_tmp $R2_tmp $R3_tmp
done


