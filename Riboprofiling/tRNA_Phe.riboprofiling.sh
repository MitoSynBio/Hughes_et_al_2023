#===================================================
# Ribo-seq script
#===================================================

# 27/06/2022


#==========================
# QUALITY CONTROL AND TRIM
# 1. TRIMGALORE
#==========================

WD=/phe_riboprofile

#### trim and FASTQC on trimmed data
TRIMGALORE=/TrimGalore-0.6.5/trim_galore
FASTQ=$WD/Data
OUT=$WD/Pipeline/trim

cd $FASTQ

for i in ribo*_R1.fastq.gz;do
	SAMPLE=${i%%_*}
	echo $SAMPLE

	$TRIMGALORE --fastqc --paired --clip_R1 4 --clip_R2 4 --three_prime_clip_R1 4 --three_prime_clip_R2 4 --gzip --output_dir $OUT ${SAMPLE}*_R1.fastq.gz ${SAMPLE}*_R2.fastq.gz
done


#==========================
# ALIGNMENT
# 1. pre-alignment
# 2. BOWTIE2 alignment
# 3. STAR alignment
#==========================

WD=/media/Data/d.rudler/phe_riboprofile

#### make bowtie index for ncRNA
REFS=/media/Data/d.rudler/species_refs/GENCODE_GRCm39.v29
cd $REFS
bowtie2-build $REFS/gencode.vM29.ncRNA_transcripts.fa ncRNA_transcripts


#### pre-align with bowtie2 against ncRNA reference genome
WD=/phe_riboprofile
INDEX=/GENCODE_GRCm39.v29/ncRNA_transcripts
INPUT=$WD/Pipeline/trim
LOG=$WD/Pipeline/pre-align/log
OUTPUT=$WD/Pipeline/pre-align

cd $INPUT

for i in ribo*R1_val_1.fq.gz;do
	SAMPLE=${i%%_*}
	echo $SAMPLE

	bowtie2 -p 10 --un-conc-gz ${OUTPUT}/${SAMPLE}_unalign -x $INDEX -1 ${SAMPLE}_*R1_val_1.fq.gz -2 $INPUT/${SAMPLE}_*R2_val_2.fq.gz 2> ${LOG}/${SAMPLE}.log | samtools sort -@ 4 -O BAM -o ${OUTPUT}/${SAMPLE}_ncRNA.bam
done


#### align with STAR
# make STAR index
cd /GENCODE_GRCm39.v29
INDEX=STAR_index_encode_141nt
FASTA=Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
GTF=Mus_musculus.GRCm39.106.gff3
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $INDEX --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang 141

# STAR fasta ref align
WD=/phe_riboprofile
INPUT=$WD/Pipeline/pre-align
OUTDIR=$WD/Pipeline/STAR_align_encode
INDEX=/GENCODE_GRCm39.v29/STAR_index_encode_141nt

cd $INPUT

for i in *unalign.1;do
	SAMPLE=${i%%_*}
	echo $SAMPLE

	STAR --runThreadN 10 --runMode alignReads --readFilesIn ${SAMPLE}_unalign.1 ${SAMPLE}_unalign.2 --readFilesCommand zcat --genomeDir $INDEX --outFileNamePrefix $OUTDIR/${SAMPLE}.STAR. --outSAMtype BAM SortedByCoordinate --alignEndsType Extend5pOfRead1 --outFilterMatchNminOverLread 0.9 --outFilterMultimapNmax 3 --alignIntronMax 2500
	samtools index -@ 8 $OUTDIR/${SAMPLE}.STAR.Aligned.sortedByCoord.out.bam
done


#==========================
# STATS
# 1. fivepseq
#==========================

WD=/phe_riboprofile
INPUT=$WD/Pipeline/STAR_align
OUTPUT=$WD/Pipeline/coverage

cd $INPUT

#### fivepseq
WD=/phe_riboprofile
BAM=$WD/Pipeline/bowtie_align
REF=/GENCODE_GRCm39.v29
OUTPUT=$WD/Pipeline/fivepseq

fivepseq -b "bed/*.bam" -g ${REF}/GRCm39.primary_assembly.genome.fa -a ${REF}/Mus_musculus.GRCm39.106.gff3 -o ${OUTPUT} -transcript-type protein_coding
