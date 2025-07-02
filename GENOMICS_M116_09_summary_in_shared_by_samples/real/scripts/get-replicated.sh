#!/bin/bash
#
#SBATCH --partition=haswell  # haswell-p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J get-rep
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-00:59# time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


IDD=M116
#CHR=chr12
CHR_LIST=($(cat ../02-vep-annot/real/inputs/chrs.list | grep -v chrY))
INPUTS_DIR=./real/inputs
OUTPUTS_DIR=./real/tmp


rm -rf ./real/tmp/*
#rm -rf ./real/tmp/M116.allchrs.var.rep

for CHR in ${CHR_LIST[@]}
do
cat ${INPUTS_DIR}/${IDD}.${CHR}.var.present | awk '{ if ($5!=1) print $1 }' | grep -v variant | grep -v "*" | grep -v chrY  > ${OUTPUTS_DIR}/${IDD}.${CHR}.var.rep
done

# use var.rep, variants replicated in at least 2 of the 3 M116´s samples, to filter results before functional analysis

cat ${OUTPUTS_DIR}/${IDD}.*.var.rep > ${OUTPUTS_DIR}/${IDD}.allchrs.var.rep

