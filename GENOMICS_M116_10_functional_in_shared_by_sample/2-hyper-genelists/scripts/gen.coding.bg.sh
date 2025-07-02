#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J coding
#SBATCH --mem 32G # memory pool for all cores
#SBATCH -t 0-03:29# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address



#cat ../../09-summary-in-shared-by-samples/real/outputs/*/M116/*.M116.*.eur_rare_all.fulltable | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM | grep Transcript  | grep CODING-GENE | grep -v -w ID | cut -f16 | sort -u | grep -v -w "" > ./real/inputs/BG_M116.coding.genes

cat ../../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/*/*.eur_rare_all.fulltable | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM | grep Transcript  | grep CODING-GENE | grep -v -w ID | cut -f16 | sort -u | grep -v -w "" > ./real/inputs/BG_1000GIBS.coding.genes


#cat ../../09-summary-in-shared-by-samples/real/outputs/*/M116/*.eur_rare_coding.genes | sort -u > ./real/inputs/M116.coding_coding.txt

#cat ../../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/*/*.eur_rare_coding.genes | sort -u > ./real/inputs/1000GIBS.coding_coding.txt


