#!/bin/bash
#
#SBATCH --partition=haswell  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J plot
#SBATCH --mem 64G # memory pool for all cores
#SBATCH -t 0-00:59# time (D-HH:MM)
#SBATCH -o ./log.%j.out # STDOUT
#SBATCH -e ./log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address




ARRAY_LIST=$1
CAT_FILE=$2
EXCL_BOOLEAN=$3
SETNAME=$4

# Array Dependent Config
# If there is no array number assigned:
#   - that means that the script is being executed directly. 
#   - it takes the first item of the list as a parameter
# If there is an array number assigned:
#   - that means that the script is being executed through the cluster
#   - it takes the item correspondent to the array number
if [ "${SLURM_ARRAY_TASK_ID}" == "" ]
then
   LEVEL=$(cat ${ARRAY_LIST} | sed -n 1p )
else
  LEVEL=$(cat ${ARRAY_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
fi 




#DB=$(ls real/outputs | grep damaging | cut -d '_' -f2 | cut -d '.' -f1)
#NUM_DB=$(ls real/outputs | grep damaging | wc -l)

# per database  a facet


#module load R/3.5.1-foss-2018b

module load R

Rscript ./real/scripts/plot.R ${LEVEL} ${CAT_FILE} ${EXCL_BOOLEAN} ${SETNAME}
