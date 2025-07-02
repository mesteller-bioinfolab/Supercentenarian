#!/bin/bash

# Main Code --------------------------------------------------------


main() {

#rm -rf real/inputs/*
#rm -rf real/outputs/*
#rm -rf real/tmp/*


TISSUE=blood


# Configuration

# CHR subset
cat ../02-vep-annot/real/inputs/chrs.list  > real/inputs/chrs.list

# All CHR
#cat ../4-vep-annot/real/inputs/chrs.list | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM > real/inputs/chrs.list

CHR_LIST=real/inputs/chrs.list


INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/06-make-table-split/real/outputs

TMP_DIR=./real/tmp

mkdir ${TMP_DIR}


IDS=$(ls ${INPUT_DIR}/${TISSUE})
#IDS=M116

mkdir ${TMP_DIR}/${TISSUE}


# eg for header
ONE_CHR_FILE=${INPUT_DIR}/${TISSUE}/M116/${TISSUE}.M116.chr22.table.split
cat ${ONE_CHR_FILE} | head -n1 > ./real/inputs/header


# Execution
# one job per ID in tissue

for ID in ${IDS}
do

echo ${ID}

mkdir ${TMP_DIR}/${TISSUE}/${ID}

  # COMMAND
  COMMAND=" \
	./real/scripts/split1000.sh \
	${CHR_LIST} \
	${INPUT_DIR} \
	${TMP_DIR} \
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




