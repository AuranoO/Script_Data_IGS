#!/bin/bash
#
#SBATCH --job-name=bbduck
#SBATCH --ntasks=1                   # nb de tâches à lancer en // (en général 1)
#SBATCH --nodes=1                    # nombre de nœuds à réserver pour chaque tâche
#SBATCH --cpus-per-task=25
#SBATCH --nodelist="jabba"

BBDUCKPATH="/home/oar-jobs/toumi/bbmap"
READPATH="/mnt/imm-storage/Data/omics/waiting_room_tmp/X204SC23021061-Z01-F001/X204SC23021061-Z01-F001/01.RawData/CedratK"
OUTPUT="/home/oar-jobs/toumi/read_trim"

${BBDUCKPATH}/bbduk.sh in1=${READPATH}/CedratK_EKDN230006071-1A_HTVMHDSX5_L4_1.fq.gz \
	 	       in2=${READPATH}/CedratK_EKDN230006071-1A_HTVMHDSX5_L4_2.fq.gz \
	 	       out1=${OUTPUT}/CedratK_trimq20_1.fq.gz out2=${OUTPUT}/CedratK_trimq20_2.fq.gz \
	 	       ref=adapter.fa \
	 	       ktrim=l \
	 	       qtrim=rl \
		       trimq=20
