#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  
  rm -rf ./real/outputs/* 

  # get samples id
  cat ./real/inputs/igsr-ibs.tsv.tsv | grep female | cut -f1 > ./real/inputs/ibs_females.txt

  SAMPLES=./real/inputs/ibs_females.txt

  TISSUE=1000GIBS

  # Configuration

  cat ./real/inputs/chrs.list | grep chrX > ./real/inputs/chrs.list1

  CHR_LIST=./real/inputs/chrs.list1
 

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/00-data/CCDG_14151_B01_GRM_WGS_2020-08-05_unnannotated
  

  #INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/00-data/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_recalibrated_variants.annotated
  # 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.annotated.vcf.gz
  TMP_DIR=./real/tmp
  OUTPUT_DIR=./real/outputs
  
  mkdir ${TMP_DIR}/splitdir
  

   #IDS=$(ls ${INPUT_DIR} | grep -v tbi)

  # Execution 
  mkdir ${OUTPUT_DIR}/${TISSUE}
  OUTPUT_DIR_THEN=${OUTPUT_DIR}/${TISSUE}

  # COMMAND
  COMMAND=" \
    ./real/scripts/create_inputs.sh \
      ${CHR_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR_THEN} \
      ${TMP_DIR} \
      ${SAMPLES}
  "


  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}


 
 exit


  # Direct execution
  #eval bash ${COMMAND}
  #exit 

  # Cluster execution
  #eval sbatch ${COMMAND}
  #exit             

} 


# FUNCTIONS ==========================================================

main
