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
	outdir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/4_peakCalling/new_analysis/${tissue}

	WT_R1=$wd/pipeline/2_alignment/${tissue}/S1_WT_${tissue}/S1_WT_${tissue}.bowtie2.nucleosome_free.175nt.bam
	WT_R2=$wd/pipeline/2_alignment/${tissue}/S3_WT_${tissue}/S3_WT_${tissue}.bowtie2.nucleosome_free.175nt.bam
	WT_R3=$wd/pipeline/2_alignment/round_2/${tissue}/22_WT4_${tissue}/22_WT4_${tissue}.bowtie2.nucleosome_free.175nt.bam
	WT_R4=$wd/pipeline/2_alignment/round_2/${tissue}/24_WT6_${tissue}/24_WT6_${tissue}.bowtie2.nucleosome_free.175nt.bam
	R1_tmp=$wd/pipeline/2_alignment/${tissue}/S1_WT_${tissue}/S1_WT_${tissue}.bowtie2.tmp.bam
	R2_tmp=$wd/pipeline/2_alignment/${tissue}/S3_WT_${tissue}/S3_WT_${tissue}.bowtie2.tmp.bam
	R3_tmp=$wd/pipeline/2_alignment/round_2/${tissue}/22_WT4_${tissue}/22_WT4_${tissue}.bowtie2.tmp.bam
	R4_tmp=$wd/pipeline/2_alignment/round_2/${tissue}/24_WT6_${tissue}/24_WT6_${tissue}.bowtie2.tmp.bam

	WT_peaks=$outdir/WT_${tissue}.forPhe11.genrich.hc.NA.nfr.narrowPeak
	WT_log=$outdir/WT_${tissue}.forPhe11.genrich.hc.NA.nfr.log
	WT_pu=$outdir/WT_${tissue}.forPhe11.genrich.pileup.hc.NA.nfr.bedgraph

	samtools sort -@ 12 -n -O BAM -o $R1_tmp $WT_R1
	samtools sort -@ 12 -n -O BAM -o $R2_tmp $WT_R2
	samtools sort -@ 12 -n -O BAM -o $R3_tmp $WT_R3
	samtools sort -@ 12 -n -O BAM -o $R4_tmp $WT_R4

	$genrich -t $R1_tmp,$R2_tmp,$R3_tmp,$R4_tmp -o $WT_peaks -E $blacklist -k $WT_pu -a 200 -m 30 -q 0.05 -j -r -e chrM,chrY -v 2> $WT_log
	gzip $WT_pu

	rm $R1_tmp $R2_tmp $R3_tmp $R4_tmp
done

## WT for Phe21
for tissue in liver
do
	#WT
	blacklist=$wd/pipeline/2_alignment/ENCFF547MET.bed
	outdir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/4_peakCalling/new_analysis/${tissue}

	WT_R1=$wd/pipeline/2_alignment/${tissue}/S1_WT_${tissue}/S1_WT_${tissue}.bowtie2.nucleosome_free.175nt.bam
	WT_R2=$wd/pipeline/2_alignment/${tissue}/S3_WT_${tissue}/S3_WT_${tissue}.bowtie2.nucleosome_free.175nt.bam
	WT_R3=$wd/pipeline/2_alignment/round_2/${tissue}/24_WT6_${tissue}/24_WT6_${tissue}.bowtie2.nucleosome_free.175nt.bam
	R1_tmp=$wd/pipeline/2_alignment/${tissue}/S1_WT_${tissue}/S1_WT_${tissue}.bowtie2.tmp.bam
	R2_tmp=$wd/pipeline/2_alignment/${tissue}/S3_WT_${tissue}/S3_WT_${tissue}.bowtie2.tmp.bam
	R3_tmp=$wd/pipeline/2_alignment/round_2/${tissue}/24_WT6_${tissue}/24_WT6_${tissue}.bowtie2.tmp.bam

	WT_peaks=$outdir/WT_${tissue}.forPhe21.genrich.hc.NA.nfr.narrowPeak
	WT_log=$outdir/WT_${tissue}.forPhe21.genrich.hc.NA.nfr.log
	WT_pu=$outdir/WT_${tissue}.forPhe21.genrich.pileup.hc.NA.nfr.bedgraph

	samtools sort -@ 12 -n -O BAM -o $R1_tmp $WT_R1
	samtools sort -@ 12 -n -O BAM -o $R2_tmp $WT_R2
	samtools sort -@ 12 -n -O BAM -o $R3_tmp $WT_R3

	$genrich -t $R1_tmp,$R2_tmp,$R3_tmp -o $WT_peaks -E $blacklist -k $WT_pu -a 200 -m 30 -q 0.05 -j -r -e chrM,chrY -v 2> $WT_log
	gzip $WT_pu

	rm $R1_tmp $R2_tmp $R3_tmp
done

for tissue in liver
do
	#Phe11
	blacklist=$wd/pipeline/2_alignment/ENCFF547MET.bed
	outdir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/4_peakCalling/new_analysis/${tissue}

	Phe11_R1=$wd/pipeline/2_alignment/${tissue}/S4_Phe11_${tissue}/S4_Phe11_${tissue}.bowtie2.nucleosome_free.175nt.bam
	Phe11_R2=$wd/pipeline/2_alignment/${tissue}/S6_Phe11_${tissue}/S6_Phe11_${tissue}.bowtie2.nucleosome_free.175nt.bam
	Phe11_R3=$wd/pipeline/2_alignment/round_2/${tissue}/20_KO5_1_1_${tissue}/20_KO5_1_1_${tissue}.bowtie2.nucleosome_free.175nt.bam
	Phe11_R4=$wd/pipeline/2_alignment/round_2/${tissue}/21_KO6_1_1_${tissue}/21_KO6_1_1_${tissue}.bowtie2.nucleosome_free.175nt.bam
	R1_tmp=$wd/pipeline/2_alignment/${tissue}/S4_Phe11_${tissue}/S4_Phe11_${tissue}.bowtie2.tmp.bam
	R2_tmp=$wd/pipeline/2_alignment/${tissue}/S6_Phe11_${tissue}/S6_Phe11_${tissue}.bowtie2.tmp.bam
	R3_tmp=$wd/pipeline/2_alignment/round_2/${tissue}/20_KO5_1_1_${tissue}/20_KO5_1_1_${tissue}.bowtie2.tmp.bam
	R4_tmp=$wd/pipeline/2_alignment/round_2/${tissue}/21_KO6_1_1_${tissue}/21_KO6_1_1_${tissue}.bowtie2.tmp.bam

	Phe11_peaks=$outdir/Phe11_${tissue}.genrich.hc.NA.nfr.narrowPeak
	Phe11_log=$outdir/Phe11_${tissue}.genrich.hc.NA.nfr.log
	Phe11_pu=$outdir/Phe11_${tissue}.genrich.pileup.hc.NA.nfr.bedgraph

	samtools sort -@ 12 -n -O BAM -o $R1_tmp $Phe11_R1
	samtools sort -@ 12 -n -O BAM -o $R2_tmp $Phe11_R2
	samtools sort -@ 12 -n -O BAM -o $R3_tmp $Phe11_R3
	samtools sort -@ 12 -n -O BAM -o $R4_tmp $Phe11_R4

	$genrich -t $R1_tmp,$R2_tmp,$R3_tmp,$R4_tmp -o $Phe11_peaks -E $blacklist -k $Phe11_pu -a 200 -m 30 -q 0.05 -j -r -e chrM,chrY -v 2> $Phe11_log
	gzip $Phe11_pu

	rm $R1_tmp $R2_tmp $R3_tmp $R4_tmp
done

for tissue in liver
do
	#Phe21
	blacklist=$wd/pipeline/2_alignment/ENCFF547MET.bed
	outdir=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq/pipeline/4_peakCalling/new_analysis/${tissue}

	Phe21_R1=$wd/pipeline/2_alignment/${tissue}/S7_Phe21_${tissue}/S7_Phe21_${tissue}.bowtie2.nucleosome_free.175nt.bam
	Phe21_R2=$wd/pipeline/2_alignment/${tissue}/S8_Phe21_${tissue}/S8_Phe21_${tissue}.bowtie2.nucleosome_free.175nt.bam
	Phe21_R3=$wd/pipeline/2_alignment/${tissue}/S9_Phe21_${tissue}/S9_Phe21_${tissue}.bowtie2.nucleosome_free.175nt.bam
	R1_tmp=$wd/pipeline/2_alignment/${tissue}/S7_Phe21_${tissue}/S7_Phe21_${tissue}.bowtie2.tmp.bam
	R2_tmp=$wd/pipeline/2_alignment/${tissue}/S8_Phe21_${tissue}/S8_Phe21_${tissue}.bowtie2.tmp.bam
	R3_tmp=$wd/pipeline/2_alignment/${tissue}/S9_Phe21_${tissue}/S9_Phe21_${tissue}.bowtie2.tmp.bam

	Phe21_peaks=$outdir/Phe21_${tissue}.genrich.hc.NA.nfr.narrowPeak
	Phe21_log=$outdir/Phe21_${tissue}.genrich.hc.NA.nfr.log
	Phe21_pu=$outdir/Phe21_${tissue}.genrich.pileup.hc.NA.nfr.bedgraph

	samtools sort -@ 12 -n -O BAM -o $R1_tmp $Phe21_R1
	samtools sort -@ 12 -n -O BAM -o $R2_tmp $Phe21_R2
	samtools sort -@ 12 -n -O BAM -o $R3_tmp $Phe21_R3

	$genrich -t $R1_tmp,$R2_tmp,$R3_tmp -o $Phe21_peaks -E $blacklist -k $Phe21_pu -a 200 -m 30 -q 0.05 -j -r -e chrM,chrY -v 2> $Phe21_log
	gzip $Phe21_pu

	rm $R1_tmp $R2_tmp $R3_tmp
done




