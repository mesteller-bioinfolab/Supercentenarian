# grep -v  -f all 1000GIBS fdr terms per category to filter results to unique ones to plot.



IBS_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/10-functional/3-webgestalt-db-gmt
CAT_FILE=./real/inputs/versions.txt
#SIG_LEVEL=nominal
#SIG_LEVEL=fdr

LEVELS=$(cat real/inputs/levels.txt  | cut -d '.' -f1)


for LEVEL in ${LEVELS[@]}

do


COMMAND="./real/scripts/excl_terms.sh ${CAT_FILE} ${IBS_DIR} ${LEVEL}"





  # Execution
  # Cluster array execution
  JOBS_COUNT=$(cat ${CAT_FILE} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  #exit





done 
