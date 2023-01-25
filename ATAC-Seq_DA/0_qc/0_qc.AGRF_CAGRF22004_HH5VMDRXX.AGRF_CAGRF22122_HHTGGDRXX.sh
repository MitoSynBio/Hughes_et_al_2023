
WD=/media/Data/s.siira/tRNA_Phe/tRNA_Phe_ATAC-Seq
DATA=$WD/data/links

## Soft links to files were created to simplify scheme naming for downstream analysis. Sample names correspond to:
# S1_WT in brain = 4_WT1_brain
# S2_WT in brain = 5_WT2_brain
# S3_WT in brain = 6_WT3_brain
# S4_Phe11 in brain = 1_KO1_1_1_brain
# S5_Phe11 in brain = 2_KO2_1_1_brain
# S6_Phe11 in brain = 3_KO4_1_1_brain
# S7_Phe21 in brain = 8_KO1_2_1_brain
# S8_Phe21 in brain = 9_KO2_2_1_brain
# S9_Phe21 in brain = 10_KO4_2_1_brain
# S10_Phe31 in brain = 1_KO1_3_1_brain
# S11_Phe31 in brain = 2_KO2_3_1_brain
# S12_Phe31 in brain = 3_KO3_3_1_brain

# S1_WT in liver = 10_WT1_liver
# S2_WT in liver = 11_WT2_liver
# S3_WT in liver = 12_WT3_liver
# S4_Phe11 in liver = 7_KO1_1_1_liver
# S5_Phe11 in liver = 8_KO2_1_1_liver
# S6_Phe11 in liver = 9_KO3_1_1_liver
# S7_Phe21 in liver = 7_KO1_2_1_liver
# S8_Phe21 in liver = 8_KO2_2_1_liver
# S9_Phe21 in liver = 9_KO3_2_1_liver
# S10_Phe31 in liver = 4_KO1_3_1_liver
# S11_Phe31 in liver = 5_KO2_3_1_liver
# S12_Phe31 in liver = 6_KO3_3_1_liver

for SAMPLE in S1_WT S2_WT S3_WT S4_Phe11 S5_Phe11 S6_Phe11 S7_Phe21 S8_Phe21 S9_Phe21 S10_Phe31 S11_Phe31 S12_Phe31 
do
	for TISSUE in brain liver
	do
		for REPLICATE in L001 L002
		do
			R1=$DATA/${TISSUE}/${SAMPLE}_${TISSUE}_${REPLICATE}_R1.fastq.gz
			R2=$DATA/${TISSUE}/${SAMPLE}_${TISSUE}_${REPLICATE}_R2.fastq.gz
			OUTDIR=$WD/pipeline/0_qc/${TISSUE}

			fastqc -t 2 -o $OUTDIR $R1 $R2
		done
	done
done


