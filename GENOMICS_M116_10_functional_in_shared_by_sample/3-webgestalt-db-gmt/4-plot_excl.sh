CATS=./real/inputs/versions.txt


LEVEL_LIST=./real/inputs/levels.txt

EXCL=TRUE

NAMESET="all"


  # COMMAND
  COMMAND=" \
    ./real/scripts/plot.sh \
    ${LEVEL_LIST} \
    ${CATS} \
    ${EXCL} \
    ${NAMESET}
  "
  
   
  # Execution 
  # Cluster array execution
  JOBS_COUNT=$(cat ${LEVEL_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
 
