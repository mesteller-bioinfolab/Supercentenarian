#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J overlong
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-03:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
module load R/3.6.0-foss-2018b

# Config

INPUTS=$1
OUTPUTS=$2
TISSUED=$3
ARRAY_LIST=$4
GENELISTS_FOLDER=$5
CAT_FILE=$6
AUTOSOMAL_GENES=$7


# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
  IDD=$(cat ${ARRAY_LIST} | sed -n 1p)
else
  IDD=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 


CATEGORIES=( $(cat ${CAT_FILE} ) )


GENELISTS=( $(ls ${GENELISTS_FOLDER} | grep txt |  cut -d '.' -f1  | rev | cut -d'_' -f2,3,4 | rev  ) )





for GENELIST in ${GENELISTS[@]}

do


for CATEGORY in ${CATEGORIES[@]}

do

# different prefix for rare cat, no additional word in filename

all="all"

if [[ "$CATEGORY" == "$all" ]]
then
        #CATEGORY_GENES=${INPUTS}/${TISSUED}.${IDD}.eur_rare.genes
    	#CATEGORY_GENES=${INPUTS}/${TISSUED}.${IDD}.eur_rare.autosomal.genes
        CATEGORY_GENES=${INPUTS}/${TISSUED}.${IDD}.eur_rare
        echo "all"
else
    	#CATEGORY_GENES=${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.genes
        #CATEGORY_GENES=${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}.autosomal.genes
        CATEGORY_GENES=${INPUTS}/${TISSUED}.${IDD}.eur_rare_${CATEGORY}
fi



cat ${CATEGORY_GENES}.genes | grep -f ${GENELISTS_FOLDER}/${GENELIST}_ensg.txt > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.${GENELIST}



#cat ${CATEGORY_GENES}.autosomal.genes | grep -f ${GENELISTS_FOLDER}/${GENELIST}_ensg.txt > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.${GENELIST}
#comm -3 ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.${GENELIST}.25chrs ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.${GENELIST} > ${OUTPUTS}/${TISSUED}.${IDD}.${CATEGORY}.${GENELIST}.chrX




done



done


