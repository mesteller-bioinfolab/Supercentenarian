#!/bin/bash

# Main Code --------------------------------------------------------


main() {

#rm -rf real/inputs/*
#rm -rf real/outputs/*
#rm -rf real/tmp/*


TISSUE=1000GIBS
INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/06-make-table-split/real/outputs
TMP_DIR=./real/tmp

# eg for header
ONE_CHR_FILE=${INPUT_DIR}/${TISSUE}/HG01510/${TISSUE}.HG01510.chr22.table.split
cat ${ONE_CHR_FILE} | head -n1 > ./real/inputs/header


#CHR_ID_LIST=../02-filter-rare/pending

CHR_ID_LIST=chr_id_list_1

# Execution
# one job per ID in tissue

  # COMMAND
  COMMAND=" \
	./real/scripts/split1000_pending.sh \
	${CHR_ID_LIST} \
	${INPUT_DIR} \
	${TMP_DIR} \
	${TISSUE}
	"


  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_ID_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}

 
exit

}

main




