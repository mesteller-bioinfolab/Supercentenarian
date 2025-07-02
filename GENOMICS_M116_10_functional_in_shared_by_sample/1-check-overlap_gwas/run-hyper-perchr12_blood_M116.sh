# Modules

module load R/3.6.0-foss-2018b

# Config


INPUTS_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/7-summary/real/outputs
OUTPUTS_DIR=./real/interactive
TISSUED=blood
IDD=M116


# Chromosomes
CHR=chr12



Rscripts ./real/scripts/hyper.R ${ESTO}
