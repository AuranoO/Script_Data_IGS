#!/bin/bash
# 
#SBATCH --job-name=RSEMBUILD
#SBATCH --ntasks=1                   # nb de tâches à lancer en // (en général 1)
#SBATCH --nodes=1                    # nombre de nœuds à réserver pour chaque tâche
#SBATCH --cpus-per-task=32           # nombre de cœurs réservés pour chaque tâche
#SBATCH --nodelist="C642002"         # check sview pour vérifier idle node

start=`date +%s`
# Préparation et Variables
GENOME="/home/oar-jobs/toumi/reference/Cedratvirus_kamtchatka_polish_bw.fasta"
INDEX="rsem_ref/ck"
GFF="/home/oar-jobs/toumi/reference/Cedratvirus_kamtchatka_filt_clean.gff"
mkdir rsem_ref

# Creating the RSEM reference using STAR
/home/oar-jobs/tools/src/RSEM/RSEM/rsem-prepare-reference --gff3 ${GFF} \
		       --gff3-genes-as-transcripts \
		       -p ${SLURM_CPUS_PER_TASK} \
		       --star \
		       --star-path /home/oar-jobs/tools/src/STAR/bin/Linux_x86_64 \
		       ${GENOME} \
	      	       ${INDEX}

# END
echo "****************************************************************"
end=`date +%s`
runtime=$(((end-start)/60))
echo ${runtime}
echo "DONE"
