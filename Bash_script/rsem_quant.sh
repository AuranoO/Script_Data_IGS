#!/bin/bash
#
#SBATCH --job-name=RSEM-QUANT
#SBATCH --ntasks=1                   # nb de tâches à lancer en // (en général 1)
#SBATCH --nodes=1                    # nombre de nœuds à réserver pour chaque tâche
#SBATCH --cpus-per-task=30           # nombre de cœurs réservés pour chaque tâche
#SBATCH --nodelist="C642001"         # check sview pour vérifier idle node

start=`date +%s`
# Variable and preparation
STARPATH="/home/oar-jobs/tools/src/STAR/bin/Linux_x86_64"
RSEMPATH="/home/oar-jobs/tools/src/RSEM/RSEM"
BAMPATH="/home/oar-jobs/toumi/mapping/pandora_merge"
INDEX="/home/oar-jobs/toumi/RSEM/index_pandora/pandora"
OUTPUT="/home/oar-jobs/toumi/RSEM/result_pandora_all"
GENOME="/home/oar-jobs/toumi/reference/A_castellanii_Neff+mitochondrion+pandora_neocaledonia.fasta"
GFF="/home/oar-jobs/toumi/reference/A_castellanii_Neff+mitochondrion+pandora_neocaledonia.gff"

mkdir -p index_pandora
mkdir -p ${OUTPUT}

# --gff3 or --gtf
${RSEMPATH}/rsem-prepare-reference --gff3 ${GFF} \
				   -p ${SLURM_CPUS_PER_TASK} \
				   --star \
				   --star-path ${STARPATH} \
				   ${GENOME} \
				   ${INDEX}

# Creating a list of my prefix for my reads
fileprefix=("A_1h" "A_3h" "A_6h" "B_1h" "B_3h" "B_6h" "C_1h" "C_3h" "C_6h" "C_1h30" "C_2h" "C_2h30" "C_4h" "C_8h" "A_MOCK" "B_MOCK" "C_MOCK")

for prefix in ${fileprefix[@]}; do
	${RSEMPATH}/rsem-calculate-expression -p ${SLURM_CPUS_PER_TASK} \
                                              --alignments \
					      --strandedness reverse \
                                              --paired-end \
					      --no-bam-output \
                                              ${BAMPATH}/${prefix}_toTranscriptome_merged.bam \
                                              ${INDEX} \
                                              ${OUTPUT}/${prefix}
        echo ${prefix}
done


echo "****************************************************************"
end=`date +%s`
runtime=$(((end-start)/60))
echo ${runtime}
echo "DONE"
