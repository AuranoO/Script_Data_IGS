#!/bin/bash
#
#SBATCH --job-name=RSEM-ALGN
#SBATCH --ntasks=1                   # nb de tâches à lancer en // (en général 1)
#SBATCH --nodes=1                    # nombre de nœuds à réserver pour chaque tâche
#SBATCH --cpus-per-task=35           # nombre de cœurs réservés pour chaque tâche
#SBATCH --nodelist="titan"         # check sview pour vérifier idle node

start=`date +%s`
# Variable and preparation
STARPATH="/home/oar-jobs/tools/src/STAR/bin/Linux_x86_64"
RSEMPATH="/home/oar-jobs/tools/src/RSEM/RSEM"
GTF="/home/oar-jobs/toumi/reference/A_castellanii_HiC-Seq+mito+haplo2m2.gtf"
READPATH="/tshared/data/rigou/Cedrat_kamchatka_transcriptome/Clean_data/"
INDEX="/home/oar-jobs/toumi/RSEM/rsem_finaltest/all"
OUTPUT="/home/oar-jobs/toumi/RSEM/result_finaltest"
GENOME="/home/oar-jobs/toumi/reference/A_castellanii_HiC-Seq+mito+haplo2m2.fasta"
mkdir -p ${OUTPUT}

# RSEM quantification
# Creating a list of my prefix for my reads
fileprefix=("A_3h" "B_3h" "C_3h" "A_MOCK" "B_MOCK" "C_MOCK")

/home/oar-jobs/tools/src/RSEM/RSEM/rsem-prepare-reference --gtf ${GTF} \
                       -p ${SLURM_CPUS_PER_TASK} \
                       --star \
                       --star-path /home/oar-jobs/tools/src/STAR/bin/Linux_x86_64 \
                       ${GENOME} \
                       ${INDEX}

for prefix in ${fileprefix[@]}; do
	${RSEMPATH}/rsem-calculate-expression -p ${SLURM_CPUS_PER_TASK} \
					      --paired-end \
					      --star \
					      --star-path ${STARPATH} \
					      --star-gzipped-read-file \
					      --paired-end \
					      --sort-bam-by-coordinate \
					      --strandedness reverse \
					      ${READPATH}/${prefix}_cleanbb_1.fq.gz \
					      ${READPATH}/${prefix}_cleanbb_2.fq.gz \
					      ${INDEX} \
					      ${OUTPUT}/${prefix}_

	mkdir ${OUTPUT}/${prefix}
	mv *${prefix}* ${prefix}/
	echo ${prefix}
done

# END
echo "****************************************************************"
end=`date +%s`
runtime=$(((end-start)/60))
echo ${runtime}
echo "DONE"
