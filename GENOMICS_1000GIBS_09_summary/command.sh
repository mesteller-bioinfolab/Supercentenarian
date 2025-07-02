#!/bin/bash

# Main Code --------------------------------------------------------

main() {
 # rm -fR ./real/outputs/*

  # this misses the non cod catts for bg
  ls ../08-filter-category-summary/real/outputs/1000GIBS/HG01501/1000GIBS.HG01501.chr12.* | grep dput | cut -d '.' -f6 | cut -d '_' -f3,4  | grep -v all > ./real/inputs/cats.txt

  # to functional inputs folders
  cat ./real/inputs/cats.txt > ../10-functional/2-hyper-genelists/real/inputs/categories.txt
  cat ./real/inputs/cats.txt > ../10-functional/2-hyper-papers//real/inputs/categories.txt



  #CATS=( damaging moderate modifier_cd modifier_nc altering ncnc_cadd15 cadd15 low_only )
  
   CATS=./real/inputs/cats.txt

   TISSUE=1000GIBS

  # Configuration
  
   INPUT_DIR=../08-filter-category-summary/real/outputs
   OUTPUT_DIR=./real/outputs  
  
   #cat real/inputs/id.list | grep -v HG01507 > real/inputs/id.list
   ls ../08-filter-category-summary/real/outputs/${TISSUE} | grep -v HG01507 > real/inputs/id.list


   ID_LIST=real/inputs/id.list1

   # Execution 

  # COMMAND
  COMMAND=" \
    ./real/scripts/summary-lists.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID_LIST} \
      ${CATS}
  "

  # Cluster array execution
  JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
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
