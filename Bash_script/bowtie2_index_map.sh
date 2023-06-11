#!/bin/bash
#
#SBATCH --job-name=BOWTIE2
#SBATCH --ntasks=1                   # nb de tâches à lancer en // (en général 1)
#SBATCH --nodes=1                    # nombre de nœuds à réserver pour chaque tâche
#SBATCH --cpus-per-task=32           # nombre de cœurs réservés pour chaque tâche
#SBATCH --nodelist="jabba"         # check sview pour vérifier idle node

start=`date +%s`
# Préparation et Variables
GENOME="/home/oar-jobs/toumi/reference/A_castellanii_HiC-Seq+mito+haplo1m2.fasta"
BOWTIEPATH="/home/oar-jobs/tools/src/bowtie2"
OUTPUT="/home/oar-jobs/toumi/mapping/bowtie_amibe_haplo1m2"
INDEX="/home/oar-jobs/toumi/mapping/index_bw_amibe_haplo1m2"
READPATH="/tshared/data/rigou/Cedrat_kamchatka_transcriptome/Clean_data"

mkdir -p ${INDEX}
mkdir -p ${OUTPUT}

${BOWTIEPATH}/bowtie2-build --threads ${SLURM_CPUS_PER_TASK} \
			    -f ${GENOME} \
			    ${INDEX}/all

# Creating a list of my prefix for my reads

#fileprefix=("A_3h" "A_1h" "A_6h" "B_1h" "B_3h" "B_6h" "C_1h" "C_3h" "C_6h" "C_1h30" "C_2h" "C_2h30" "C_4h" "C_8h" "A_MOCK" "B_MOCK" "C_MOCK")
fileprefix=("A_3h" "B_3h" "C_3h")

# Loop through my list and map each pair of reads
for prefix in ${fileprefix[@]}; do
	${BOWTIEPATH}/bowtie2 -x ${INDEX}/all \
			      -p ${SLURM_CPUS_PER_TASK} \
			      -1 ${READPATH}/${prefix}_cleanbb_1.fq.gz \
			      -2 ${READPATH}/${prefix}_cleanbb_2.fq.gz \
			      2> ${OUTPUT}/${prefix}_stat.file \
			      | samtools view -@ ${SLURM_CPUS_PER_TASK} -h -b - \
			      | samtools sort -o ${OUTPUT}/${prefix}_sorted.bam

	echo $prefix

done

echo "********************************************************************"
end=`date +%s`
runtime=$(((end-start)/60))
echo $runtime
echo "DONE"
