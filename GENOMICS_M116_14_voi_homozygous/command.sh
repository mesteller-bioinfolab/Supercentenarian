#!/bin/bash

# Main Code --------------------------------------------------------


main() {
   rm -rf ./real/inputs/*
   rm -rf ./real/tmp/*
   rm -rf ./real/outputs/*



   # Configuration

    TMP_DIR=./real/tmp
   

    TISSUE=2TISSUES
 
    #cat ../09-summary-in-shared-by-samples/real/outputs/*/M116/*_novel.variants > ./real/inputs/novel.txt
    #cat ../../1000GIBS/08-filter-category-summary/real/inputs/ids_already.txt > ./real/inputs/ids.txt
   

    echo M116 >> ./real/inputs/ids.txt

  
    IDS_LIST=./real/inputs/ids.txt
   
    OUTPUT_DIR=./real/outputs

   # HEADERS
   cat ../09-summary-in-shared-by-samples/real/outputs/urine/M116/urine.M116.chrM.eur_rare_all.fulltable | head -n1 > ${TMP_DIR}/2TISSUES.M116.header
   #cat ../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/HG01501/1000GIBS.HG01501.chr22.eur_rare_all.fulltable | head -n1 > ${TMP_DIR}/1000GIBS.header


 

  # COMMAND
  COMMAND=" \
    ./real/scripts/voi_homozygous.sh \
      ${TMP_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE}
      ${IDS_LIST}
  "


  # Cluster array execution
  JOBS_COUNT=$(cat ${IDS_LIST} | wc -l)
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

