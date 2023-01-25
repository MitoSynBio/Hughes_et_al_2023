wd=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq


#featureCounts v2.0.0
### This analysis used
### NFR region ####
# count reads per tRNA + flank (50nt either side) - round 1
irds=/media/IRDS/Server/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
tRNAs=$wd/pipeline/7_tRNA_specific/gtRNAdb_mm10.50nt_flank.gtf
input=`find $irds/pipeline/2_alignment/brain $irds/pipeline/2_alignment/liver -name "*.nucleosome_free.175nt.bam"`
output=$wd/pipeline/7_tRNA_specific/tRNAPhe_ATACSeq.tRNA_50ntflanks.round_1.nfr.counts
log=$wd/pipeline/7_tRNA_specific/tRNAPhe_ATACSeq.tRNA_50ntflanks.featureCounts.round_1.nfr.log

featureCounts -p -B -C -P -T 12 -d 20 -D 1000 -Q 30 -f -t 'exon' -g 'transcript_id' -a $tRNAs -o $output $input &> $log


# count reads per tRNA + flank - round 2
irds=/media/IRDS/Server/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
tRNAs=$wd/pipeline/7_tRNA_specific/gtRNAdb_mm10.50nt_flank.gtf
input=`find $irds/pipeline/2_alignment/round_2/brain $irds/pipeline/2_alignment/round_2/liver -name "*.nucleosome_free.175nt.bam"`
output=$wd/pipeline/7_tRNA_specific/tRNAPhe_ATACSeq.tRNA_50ntflanks.round_2.nfr.counts
log=$wd/pipeline/7_tRNA_specific/tRNAPhe_ATACSeq.tRNA_50ntflanks.featureCounts.round_2.nfr.log

featureCounts -p -B -C -P -T 12 -d 20 -D 1000 -Q 30 -f -t 'exon' -g 'transcript_id' -a $tRNAs -o $output $input &> $log


# count reads per tRNA + flank - round 3
irds=/media/IRDS/Server/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
tRNAs=$wd/pipeline/7_tRNA_specific/gtRNAdb_mm10.50nt_flank.gtf
input=`find $wd/pipeline/2_alignment/AGRF_CAGRF21118859_HVGVVDRXY/liver -name "*.nucleosome_free.175nt.bam"`
output=$wd/pipeline/7_tRNA_specific/tRNAPhe_ATACSeq.tRNA_50ntflanks.round_3.nfr.counts
log=$wd/pipeline/7_tRNA_specific/tRNAPhe_ATACSeq.tRNA_50ntflanks.featureCounts.round_3.nfr.log

featureCounts -p -B -C -P -T 12 -d 20 -D 1000 -Q 30 -f -t 'exon' -g 'transcript_id' -a $tRNAs -o $output $input &> $log


