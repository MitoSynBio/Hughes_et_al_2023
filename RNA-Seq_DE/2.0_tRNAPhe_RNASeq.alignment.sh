WD=/media/Data/s.siira/tRNA_Phe_RNAseq
REFDIR=/media/Data/ref_seqs/Mus_musculus/mm10_GRCm38/GENCODE
MASKED=$REFDIR/GRCm38.primary_assembly.genome.numts_masked.fa
INDEX=$REFDIR/vM24/STAR.GENCODE_vM24_w_custom_chrM.150nt.NUMTs_masked
GTF=/media/Data/ref_seqs/Mus_musculus/mm10_GRCm38/GENCODE/vM24/gencode.vM24.primary_assembly.annotation.w_custom_chrM.gtf

#mkdir $INDEX 
#STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $INDEX --genomeFastaFiles $MASKED --sjdbGTFfile $GTF --sjdbOverhang 149

for SAMPLE in 1_H2G73DSXY_CCGCGGTT-AGCGCTAG 2_H2G73DSXY_TTATAACC-GATATCGA 3_H2G73DSXY_GGACTTGG-CGCAGACG 4_H2G73DSXY_AAGTCCAA-TATGAGTA 5_H2G73DSXY_ATCCACTG-AGGTGCGT 6_H2G73DSXY_GCTTGTCA-GAACATAC 7_H2G73DSXY_CAAGCTAG-ACATAGCG 8_H2G73DSXY_TGGATCGA-GTGCGATA 9_H2G73DSXY_AGTTCAGG-CCAACAGA 10_H2G73DSXY_GACCTGAA-TTGGTGAG 11_H2G73DSXY_TCTCTACT-CGCGGTTC 12_H2G73DSXY_CTCTCGTC-TATAACCT 13_H2G73DSXY_CCAAGTCT-AAGGATGA 14_H2G73DSXY_TTGGACTC-GGAAGCAG 15_H2G73DSXY_GGCTTAAG-TCGTGACC 16_H2G73DSXY_AATCCGGA-CTACAGTT 17_H2G73DSXY_TAATACAG-ATATTCAC 18_H2G73DSXY_CGGCGTGA-GCGCCTGT 19_H2G73DSXY_ATGTAAGT-ACTCTATG 20_H2G73DSXY_GCACGGAC-GTCTCGCA 21_H2G73DSXY_GGTACCTT-AAGACGTC 22_H2G73DSXY_AACGTTCC-GGAGTACT 23_H2G73DSXY_GCAGAATT-ACCGGCCA 24_H2G73DSXY_ATGAGGCC-GTTAATTG
do
	DATA=$WD/pipeline/1_trimming/${SAMPLE}
	L1R1=$DATA/${SAMPLE}_L002_R1_val_1.fq.gz
	L1R2=$DATA/${SAMPLE}_L002_R2_val_2.fq.gz
	L2R1=$DATA/${SAMPLE}_L003_R1_val_1.fq.gz
	L2R2=$DATA/${SAMPLE}_L003_R2_val_2.fq.gz
	OUTDIR=$WD/pipeline/2_alignment/${SAMPLE}

	mkdir $OUTDIR
	STAR --runThreadN 16 --runMode alignReads --readFilesIn $L1R1,$L2R1 $L1R2,$L2R2 --readFilesCommand zcat --genomeDir $INDEX --outFileNamePrefix $OUTDIR/${SAMPLE}. --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts
	samtools index -@ 16 $OUTDIR/${SAMPLE}.Aligned.sortedByCoord.out.bam

	PICARD=/media/Data/tools/picard-tools/picard-2.22.1/picard.jar
	INPUT=$OUTDIR/${SAMPLE}.Aligned.toTranscriptome.out.bam
	METRICS=$OUTDIR/${SAMPLE}.STAR_transcriptome.metrics
	HIST=$OUTDIR/${SAMPLE}.STAR_transcriptome.insert_size_histogram.pdf

	# reads mapped in proper primary pair
	samtools view -@ 16 -b -f 2 -F 260 $INPUT | java -jar $PICARD CollectInsertSizeMetrics I=/dev/stdin O=$METRICS H=$HIST
done


TXFASTA=$WD/pipeline/2_alignment/Tx.fa
TXGFF=$WD/pipeline/2_alignment/Tx.gff
gffread $GTF -g $MASKED -w $TXFASTA -o $TXGFF


### alignment-based quantification
## Quantify
for SAMPLE in 1_H2G73DSXY_CCGCGGTT-AGCGCTAG 2_H2G73DSXY_TTATAACC-GATATCGA 3_H2G73DSXY_GGACTTGG-CGCAGACG 4_H2G73DSXY_AAGTCCAA-TATGAGTA 5_H2G73DSXY_ATCCACTG-AGGTGCGT 6_H2G73DSXY_GCTTGTCA-GAACATAC 7_H2G73DSXY_CAAGCTAG-ACATAGCG 8_H2G73DSXY_TGGATCGA-GTGCGATA 9_H2G73DSXY_AGTTCAGG-CCAACAGA 10_H2G73DSXY_GACCTGAA-TTGGTGAG 11_H2G73DSXY_TCTCTACT-CGCGGTTC 12_H2G73DSXY_CTCTCGTC-TATAACCT 13_H2G73DSXY_CCAAGTCT-AAGGATGA 14_H2G73DSXY_TTGGACTC-GGAAGCAG 15_H2G73DSXY_GGCTTAAG-TCGTGACC 16_H2G73DSXY_AATCCGGA-CTACAGTT 17_H2G73DSXY_TAATACAG-ATATTCAC 18_H2G73DSXY_CGGCGTGA-GCGCCTGT 19_H2G73DSXY_ATGTAAGT-ACTCTATG 20_H2G73DSXY_GCACGGAC-GTCTCGCA 21_H2G73DSXY_GGTACCTT-AAGACGTC 22_H2G73DSXY_AACGTTCC-GGAGTACT 23_H2G73DSXY_GCAGAATT-ACCGGCCA 24_H2G73DSXY_ATGAGGCC-GTTAATTG
do
	OUTDIR=$WD/pipeline/2_alignment/${SAMPLE}/quantification
	INPUT=$WD/pipeline/2_alignment/${SAMPLE}/${SAMPLE}.Aligned.toTranscriptome.out.bam

	mkdir $OUTDIR
	salmon quant -p 16 -t $TXFASTA -l ISR -a $INPUT --seqBias --gcBias -o $OUTDIR 
done

for SAMPLE in 1_H2G73DSXY_CCGCGGTT-AGCGCTAG 2_H2G73DSXY_TTATAACC-GATATCGA 3_H2G73DSXY_GGACTTGG-CGCAGACG 4_H2G73DSXY_AAGTCCAA-TATGAGTA 5_H2G73DSXY_ATCCACTG-AGGTGCGT 6_H2G73DSXY_GCTTGTCA-GAACATAC 7_H2G73DSXY_CAAGCTAG-ACATAGCG 8_H2G73DSXY_TGGATCGA-GTGCGATA 9_H2G73DSXY_AGTTCAGG-CCAACAGA 10_H2G73DSXY_GACCTGAA-TTGGTGAG 11_H2G73DSXY_TCTCTACT-CGCGGTTC 12_H2G73DSXY_CTCTCGTC-TATAACCT 13_H2G73DSXY_CCAAGTCT-AAGGATGA 14_H2G73DSXY_TTGGACTC-GGAAGCAG 15_H2G73DSXY_GGCTTAAG-TCGTGACC 16_H2G73DSXY_AATCCGGA-CTACAGTT 17_H2G73DSXY_TAATACAG-ATATTCAC 18_H2G73DSXY_CGGCGTGA-GCGCCTGT 19_H2G73DSXY_ATGTAAGT-ACTCTATG 20_H2G73DSXY_GCACGGAC-GTCTCGCA 21_H2G73DSXY_GGTACCTT-AAGACGTC 22_H2G73DSXY_AACGTTCC-GGAGTACT 23_H2G73DSXY_GCAGAATT-ACCGGCCA 24_H2G73DSXY_ATGAGGCC-GTTAATTG
do
	samtools idxstat -@ 12 $WD/pipeline/2_alignment/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.bam > $WD/pipeline/2_alignment/${SAMPLE}/${SAMPLE}.Aligned.sortedByCoord.out.idxstat
done

