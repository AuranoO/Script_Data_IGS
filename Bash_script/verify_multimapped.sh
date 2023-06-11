#!/bin/bash
# 
#SBATCH --job-name=multimapped-check
#SBATCH --ntasks=1                   # nb de tâches à lancer en // (en général 1)
#SBATCH --nodes=1                    # nombre de nœuds à réserver pour chaque tâche
#SBATCH --cpus-per-task=30           # nombre de cœurs réservés pour chaque tâche
#SBATCH --nodelist="C642001"         # check sview pour vérifier idle node

start=`date +%s`

fileprefix=("A_1h" "A_3h" "A_6h" "A_MOCK" "B_1h" "B_3h" "B_6h" "B_MOCK" "C_1h" "C_3h" "C_6h" "C_MOCK" "C_1h30" "C_2h" "C_2h30" "C_4h" "C_8h")
for prefix in ${fileprefix[@]}; do

	samtools view -@ ${SLURM_CPUS_PER_TASK} STAR_${prefix}Aligned.sorted* -o ${prefix}.bam.sam

	sed -i 's/scaffold_[0-9]*/amibe/g' ${prefix}.bam.sam

	#gawk 't_ass[$1 $3]++==0{t_ass2[$1]=sprintf("%s %s",t_ass2[$1],$3)}END{for(k in t_ass2){print k,t_ass2[k]}}' ${prefix}.bam.sam|gawk 'NF>2'|grep -v "*"; done

	gawk 't_ass[$1 $3]++==0{t_ass2[$1]=sprintf("%s %s",t_ass2[$1],$3)}END{for(k in t_ass2){print k,t_ass2[k]}}' ${prefix}.bam.sam|gawk 'NF>2' > ${prefix}.bam.sam.csv

	#CLEAN UP
	rm ${prefix}.bam.sam
done
# END
echo "****************************************************************"
end=`date +%s`
runtime=$(((end-start)/60))
echo ${runtime}
echo "DONE"
