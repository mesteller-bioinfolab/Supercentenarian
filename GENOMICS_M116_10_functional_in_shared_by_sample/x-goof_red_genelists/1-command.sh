#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -fR ./real/inputs/*
  rm -fR ./real/tmp/*
  rm -fR ./real/outputs/*

  # Configuration

  ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > ./real/inputs/cats.txt  
  
  # All CHRs
  cat ../../02-vep-annot/real/inputs/chrs.list | grep -v chrM | grep -v chrX | grep -v chrY > ./real/inputs/chrs.list
  CHR_LIST=./real/inputs/chrs.list
   

    INPUT_DIR=../../09-summary-in-shared-by-samples/real/outputs
    

    #INPUT_DIR=../../09-summary-in-shared-by-samples/real/outputs
    #INPUT_DIR2=../../1000GIBS/08-filter-category-summary/real/outputs
    #INPUT_DIR=../../1000GIBS/08-filter-category-summary/real/outputs
   
    
    OUTPUT_DIR=./real/inputs 

    #echo M116 > ./real/inputs/ids_list.txt
    #ls -l ../../1000GIBS/09-summary/real/outputs/ | grep txt | grep -v -w 116 | cut -d ' ' -f14 | cut -d '.' -f2 >> ./real/inputs/ids_list.txt
    #IDS=$(cat ./real/inputs/ids_list.txt)
 
  
    ID=M116
    
    #cp ../../11-g-matrix/real/inputs/cats.txt real/inputs
    CATS=./real/inputs/cats.txt
  # Execution 
  

  # COMMAND
  COMMAND=" \
    ./real/scripts/inputs.sh \
      ${CHR_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${ID} \
      ${CATS}
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
