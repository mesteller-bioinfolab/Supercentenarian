#!/bin/bash

# Main Code --------------------------------------------------------

main() {


TISSUE=1000GIBS

# Configuration

#cat ../03-vep-annot/real/inputs/chrs.list | grep -v -w chrM | grep -v -w chrX | grep -v -w chrY > ./real/inputs/chrs.list
#CHR_LIST=./real/inputs/chrs.list



INPUT_DIR=./real/tmp
OUTPUT_DIR=./real/outputs


ONE_TMP_FILE=real/tmp/1000GIBS/HG02232/1000GIBS.HG02232.chr9_aa.table.split2.insplits
cat ${ONE_TMP_FILE} | head -n1  > real/inputs/splitheader

CHR_ID_LIST=chr_id_list

#IDS=$(ls ${INPUT_DIR}/${TISSUE})

# Execution 

mkdir ${OUTPUT_DIR}/${TISSUE}

# one job per ID in tissue

#for ID in ${IDS}
#do
#echo $ID
#mkdir ${OUTPUT_DIR}/${TISSUE}/${ID}


  # COMMAND
  COMMAND=" \
	./real/scripts/rejoin1000_failed.sh \
	${CHR_ID_LIST} \
	${INPUT_DIR} \
	${OUTPUT_DIR} \
	${TISSUE}
	"


  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_ID_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}

#done
 
exit


}

main


