CATS=./real/inputs/versions.txt

LEVEL_LIST=./real/inputs/levels.txt



  # COMMAND
  COMMAND=" \
    ./real/scripts/plot.sh \
    ${LEVEL_LIST} \
    ${CATS}
  "
  
   
  # Execution 
  # Cluster array execution
  JOBS_COUNT=$(cat ${LEVEL_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
 
