
WD=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
DATA=$WD/data/links
TG=/media/Data/tools/TrimGalore/TrimGalore-0.6.5/trim_galore


for SAMPLE in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31
do
	for TISSUE in brain liver
	do
		for REPLICATE in L001 L002
		do
			R1=$DATA/${TISSUE}/${SAMPLE}_${TISSUE}_${REPLICATE}_R1.fastq.gz
			R2=$DATA/${TISSUE}/${SAMPLE}_${TISSUE}_${REPLICATE}_R2.fastq.gz
			OUTDIR=$WD/pipeline/1_trimming/${TISSUE}/${SAMPLE}_${REPLICATE}

			mkdir $OUTDIR
			$TG --paired --fastqc --output_dir $OUTDIR $R1 $R2
		done
	done
done

