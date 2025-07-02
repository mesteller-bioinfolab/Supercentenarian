#!/bin/bash
#
#SBATCH --partition=haswell  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J chrsjoin
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-02:59# time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules

module load R


# Config
ARRAY_LIST=$1
INPUTS=$2   
OUTPUTS=$3
IDD=$4


# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  CHR=$(cat ${ARRAY_LIST} | sed -n 1p)
else
  CHR=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 





# M116
cat $(ls ${INPUTS}/*/${IDD}/* | grep --color=auto genevar | grep all | grep -w ${CHR} ) | head -n1 > ./real/inputs/${IDD}.${CHR}.genevar
cat $(ls ${INPUTS}/*/${IDD}/* | grep --color=auto genevar | grep all | grep -w ${CHR} ) | grep -v -w ENSG   | sort -u  >> ./real/inputs/${IDD}.${CHR}.genevar
cat ./real/inputs/${IDD}.${CHR}.genevar | cut -f1 | sort -u | grep -v -w ENSG  > ./real/inputs/${IDD}.${CHR}.genes


# All genes with rare variants/VOIs in IBS plus M116
# with any rare variants - filter by actual EUR_unrel AF in 1000G 30X VCFs (1000G 00-data common_snvs)
cat $(ls ../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/*/*.genevar | grep -w ${CHR} |  grep all ) | cut -f1  |  grep -v -w ENSG  | sort -u > ./real/inputs/1000GIBS.${CHR}.genes
cat ./real/inputs/${IDD}.${CHR}.genes ./real/inputs/1000GIBS.${CHR}.genes | sort -u > ./real/inputs/everyones.${CHR}.genes

