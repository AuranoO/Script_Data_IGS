#!/bin/bash
#
#SBATCH --job-name=RSEM-ALGN
#SBATCH --ntasks=1                   # nb de tâches à lancer en // (en général 1)
#SBATCH --nodes=1                    # nombre de nœuds à réserver pour chaque tâche
#SBATCH --cpus-per-task=32           # nombre de cœurs réservés pour chaque tâche
#SBATCH --nodelist="C642002"         # check sview pour vérifier idle node

start=`date +%s`
# Variable and preparation
BOWTIEPATH="/home/oar-jobs/tools/src/bowtie2"
RSEMPATH="/home/oar-jobs/tools/src/RSEM/RSEM"
READPATH="/tshared/data/rigou/Cedrat_kamchatka_transcriptome/Clean_data/"
OUTPUT="/home/oar-jobs/toumi/RSEM/result_test_HiC_bowtie"
GENOME="/home/oar-jobs/toumi/reference/Neff_assembly.fa"
INDEX="index_bowtie_neff_HiC/amibe"
GTF="/home/oar-jobs/toumi/reference/Neff_annotations.gtf"
mkdir -p index_bowtie_neff_HiC
mkdir -p ${OUTPUT}

# Creating the RSEM reference using BOWTIE2
${RSEMPATH}/rsem-prepare-reference --gtf ${GTF} \
                        	   -p ${SLURM_CPUS_PER_TASK} \
                        	   --bowtie2 \
				   --polyA \
                        	   --bowtie2-path ${BOWTIEPATH} \
                        	   ${GENOME} \
                        	   ${INDEX}

# RSEM quantification
# Creating a list of my prefix for my reads
#fileprefix=("A_1h" "A_3h" "A_6h" "B_1h" "B_3h" "B_6h" "C_1h" "C_3h" "C_6h" "C_1h30" "C_2h" "C_2h30" "C_4h" "C_8h" "A_MOCK" "B_MOCK" "C_MOCK")
fileprefix=("C_MOCK")

for prefix in ${fileprefix[@]}; do
	${RSEMPATH}/rsem-calculate-expression --paired-end \
					      -p ${SLURM_CPUS_PER_TASK} \
					      --bowtie2 \
					      --bowtie2-path ${BOWTIEPATH} \
#					      --sort-bam-by-coordinate \
					      --strandedness reverse \
					      ${READPATH}/${prefix}_cleanbb_1.fq.gz \
					      ${READPATH}/${prefix}_cleanbb_2.fq.gz \
					      ${INDEX} \
					      ${OUTPUT}/${prefix}_
	echo ${prefix}
done

# END
echo "****************************************************************"
end=`date +%s`
runtime=$(((end-start)/60))
echo ${runtime}
echo "DONE"
