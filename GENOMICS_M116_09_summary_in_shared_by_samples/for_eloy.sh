#!/bin/bash
#
#SBATCH --partition=normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J summary
#SBATCH --mem 120G # memory pool for all cores
#SBATCH -t 0-19:59# time (D-HH:MM)



CHR=$1



# grepping in the join fulltables for the 2tissues variants

cat real/outputs/blood/M116/blood.M116.chr22.eur_rare_all.fulltable  | head -n1 > 2TISSUES.${CHR}.fulltable

cat real/outputs/*/M116/*.${CHR}.eur_rare_all.fulltable | sort -u | grep -w -f real/tmp/M116.${CHR}.var.rep  >> 2TISSUES.${CHR}.fulltable
