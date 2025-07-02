#!/bin/bash

# Main Code --------------------------------------------------------

main() {

TISSUE=blood

# Configuration


# Subset
#cat ../4-vep-annot/real/inputs/chrs.list | grep -w chr2 > real/inputs/chrs.list
#cat ../4-vep-annot/real/inputs/chrs.list | grep chr4 >> real/inputs/chrs.list
#cat ../4-vep-annot/real/inputs/chrs.list | grep chr5 >> real/inputs/chrs.list
#cat ../4-vep-annot/real/inputs/chrs.list | grep chr20 >> real/inputs/chrs.list
#CHR_LIST=./real/inputs/chrs.list


# Whole set
#cat ../4-vep-annot/real/inputs/chrs.list | grep -v -w chrM | grep -v -w chrX | grep -v -w chrY > ./real/inputs/chrs.list
CHR_LIST=./real/inputs/chrs.list


#INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/7-make-table-split/real/outputs

INPUT_DIR=./real/tmp
OUTPUT_DIR=./real/outputs


ONE_TMP_FILE=real/tmp/blood/M116/blood.M116.chr1_aa.table.split2.insplits
cat ${ONE_TMP_FILE} | head -n1  > real/inputs/splitheader


IDS=$(ls ${INPUT_DIR}/${TISSUE})

# Execution 

mkdir ${OUTPUT_DIR}/${TISSUE}

# one job per ID in tissue

for ID in ${IDS}

do

echo $ID

mkdir ${OUTPUT_DIR}/${TISSUE}/${ID}


  # COMMAND
  COMMAND=" \
	./real/scripts/rejoin1000.sh \
	${CHR_LIST} \
	${INPUT_DIR} \
	${OUTPUT_DIR} \
	${TISSUE} \
	${ID}
	"


  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}

done
 
exit


}

main


