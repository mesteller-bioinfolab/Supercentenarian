#!/bin/bash
#
#SBATCH --partition=haswell # haswell  #normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J cadd-split
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 9-23:59# time (D-HH:MM)
#SBATCH -o ./real/tmp/log.%j.out # STDOUT
#SBATCH -e ./real/tmp/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


IDD=$1
TISSUED=$2
CHR=$3
KEY=$4


module load Python/3.7.2-GCCcore-8.2.0


INPUTS_DIR=./real/tmp
OUTPUTS_DIR=./real/tmp

TABLE=${INPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.header


python3 ./real/scripts/split-VEP-field-CADD_2.py -i ${TABLE} -o ${OUTPUTS_DIR}/${TISSUED}/${IDD}/${TISSUED}.${IDD}.${CHR}.table.split2.insplits -k ${KEY}

