#!/bin/bash
#
#SBATCH --job-name=STARWFL
#SBATCH --ntasks=1                   # nb de tâches à lancer en // (en général 1)
#SBATCH --nodes=1                    # nombre de nœuds à réserver pour chaque tâche
#SBATCH --cpus-per-task=20           # nombre de cœurs réservés pour chaque tâche
#SBATCH --nodelist="C642002"

start=`date +%s`
# Préparation et Variables
GENOME="/home/oar-jobs/toumi/reference/A_castellanii_HiC-Seq+mito+haplo2m2.fasta"
OUTPUT="/home/oar-jobs/toumi/mapping/star_hic_mito_h2m22"
INDEX="/home/oar-jobs/toumi/mapping/index_hic_mito_h2m22"
GTF="/home/oar-jobs/toumi/reference/A_castellanii_HiC-Seq+mito+haplo2m2.gtf"
STARPATH="/home/oar-jobs/tools/src/STAR/bin/Linux_x86_64"
READPATH="/tshared/data/rigou/Cedrat_kamchatka_transcriptome/Clean_data"
# Creating the STAR index using reference genome

rm -r ${INDEX}
rm -r ${OUTPUT}
mkdir ${INDEX}
mkdir ${OUTPUT}

${STARPATH}/STAR --runMode genomeGenerate \
		 --runThreadN ${SLURM_CPUS_PER_TASK} \
		 --genomeDir ${INDEX} \
		 --genomeFastaFiles $GENOME \
		 --sjdbGTFfile $GTF \
		 --genomeSAindexNbases 11

# We can add a gtf or gff file using --sjdbGTFtagExonParentTranscript some.gff or --sjdbGTFfile some.gtf

# Creating a list of my prefix for my reads
#fileprefix=("A_1h" "A_3h" "A_6h" "A_MOCK" "B_1h" "B_3h" "B_6h" "B_MOCK" "C_1h" "C_3h" "C_6h" "C_MOCK" "C_1h30" "C_2h" "C_2h30" "C_4h" "C_8h")
#fileprefix=("T1" "T2" "T3" "T4" "T5" "T6" "T7" "T8" "T9" "T10" "T11" "T12" "T13" "T14" "T15")
fileprefix=("A_3h" "B_3h" "C_3h")

# Loop through my list and map each pair of reads
for prefix in ${fileprefix[@]}; do
	${STARPATH}/STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
			 --genomeDir ${INDEX} \
			 --readFilesCommand zcat -c \
			 --readFilesIn ${READPATH}/${prefix}_cleanbb_1.fq.gz ${READPATH}/${prefix}_cleanbb_2.fq.gz \
			 --outFileNamePrefix ${OUTPUT}/STAR_${prefix}. \
			 --outSAMtype BAM SortedByCoordinate \
			 --quantMode GeneCounts TranscriptomeSAM \
			 --twopassMode Basic \
			 --alignIntronMax 3000 \
			 --alignMatesGapMax 3000 \
#			 --outSAMunmapped Within \
#			 --limitBAMsortRAM 1006353092

	echo $prefix
samtools index ${OUTPUT}/STAR_${prefix}.Aligned.sortedByCoord.out.bam
samtools flagstat ${OUTPUT}/STAR_$prefix.Aligned.sortedByCoord.out.bam > ${OUTPUT}/STAR_$prefix.Aligned.sortedByCoord.out.flagstat
#samtools flagstat ${OUTPUT}/STAR_$prefix.Aligned.toTranscriptome.out.bam > ${OUTPUT}/STAR_$prefix.Aligned.toTranscriptome.out.flagstat
done

echo "Clean up"
cd ${OUTPUT}
rm -r *SJ* *_STARgenome *_STARpass1 *_STARtmp *Log.out *Log.progress.out

echo "****************************************************************"
end=`date +%s`
runtime=$(((end-start)/60))
echo ${runtime}
echo "DONE"
