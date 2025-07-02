CAT_FILE=./real/inputs/versions.txt

SIG_LEVEL=nominal

SIG_LEVELS=./real/inputs/levels

COMMAND="./real/scripts/poolterms.sh ${CAT_FILE} ${SIG_LEVELS}"


  # Execution
  # Cluster array execution
  JOBS_COUNT=$(cat ${CAT_FILE} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  exit




#cat $(ls ./real/outputs/*.fdr.txt | grep "/${CATEGORY}_" )  | cut -f1,2 | sort -u  > 1000GIBS.${CATEGORY}.fdr


#done

