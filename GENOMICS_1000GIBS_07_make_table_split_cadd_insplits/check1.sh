#!/bin/bash
#
#SBATCH -p haswell # normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J check1
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:59# time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address

IDS=$(ls real/tmp/1000GIBS)

echo ID list headers  > insplits_chr_step1.txt
#echo ID headers insplits > insplits_chr_step2.txt

for IDD in ${IDS}
do
for CHR	 in $(seq 1 22)
do

# step 1
echo ${IDD} ${CHR} $(cat real/inputs/1000GIBS.${IDD}.chrs_split.list  | grep chr${CHR}_ | wc -l) $(la real/tmp/1000GIBS/${IDD}/* | grep chr${CHR}_ | grep header | wc -l) >> insplits_chr_step1.txt
 
# step 2
# echo ${IDD} ${CHR} $(la real/tmp/1000GIBS/${IDD}/* | grep chr${CHR}_ | grep "insplits" | wc -l)  $(la real/tmp/1000GIBS/${IDD}/* | grep chr${CHR}_ | grep header | wc -l)>> insplits_chr_step2.txt

done
done
