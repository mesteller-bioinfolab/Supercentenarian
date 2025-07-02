#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J getlist
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-06:59# time (D-HH:MM)
#SBATCH -o ./real/inputs/log.%j.out # STDOUT
#SBATCH -e ./real/inputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



# Modules
#module load VEP/102.0-foss-2019a-Perl-5.28.1
module load BCFtools/1.16-GCC-11.2.0


# Config

INPUTS=$1
OUTPUTS=$2
ARRAY_LIST=$3



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



# get to exclude

cat  ${INPUTS}/1000G.${CHR}.eur.common.table |  cut -d ' ' -f1,2,3,4 |  tr ' ' ':' | awk '{ gsub("chr", "") ; print $0 }' > ${OUTPUTS}/1000G.${CHR}.eur.common.id
 



# filter -e # other script to not confude INPUTS DIRS
#bcftools view -e 'ID=@test.snp' real/outputs/1000GIBS/HG02222/1000GIBS.HG02222.chr1.filtered.gz
