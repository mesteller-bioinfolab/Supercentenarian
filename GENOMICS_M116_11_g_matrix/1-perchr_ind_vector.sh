#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -fR ./real/tmp/*
  
  cp ../09-summary-in-shared-by-samples/real/inputs/cats.txt real/inputs
  CATS=./real/inputs/cats.txt
  

  # Configuration

  # All CHRs
  cat ../02-vep-annot/real/inputs/chrs.list | grep -v chrM | grep -v chrX | grep -v chrY > ./real/inputs/chrs.list
  CHR_LIST=./real/inputs/chrs.list
   

   
    INPUT_DIR1=../09-summary-in-shared-by-samples/real/outputs
    INPUT_DIR2=../../1000GIBS/08-filter-category-summary/real/outputs
  
    #INPUT_DIR=../../1000GIBS/08-filter-category-summary/real/outputs
    OUTPUT_DIR=./real/tmp
  
 
    echo M116 > ./real/inputs/ids_list.txt
    #ls -l ../../1000GIBS/09-summary/real/outputs/ | grep txt | grep -v -w 116 | cut -d ' ' -f14 | cut -d '.' -f2 >> ./real/inputs/ids_list.txt
    ls ../../1000GIBS/09-summary/real/outputs/*txt | cut -d '.' -f6 >> ./real/inputs/ids_list.txt
    #IDS=$(cat ./real/inputs/ids_list.txt)
 

    ID_LIST=./real/inputs/ids_list.txt
  # Execution 
   
  # one job per ID in tissue
#for ID in ${IDS[@]}
#do
#echo ${ID}
#if [ "${ID}" == "M116" ]
#then 
#    INPUT_DIR=${INPUT_DIR1}
#else
#    INPUT_DIR=${INPUT_DIR2}
#fi




  # COMMAND
  COMMAND=" \
    ./real/scripts/voi_g_matrix.sh \
      ${ID_LIST} \
      ${INPUT_DIR1} \
      ${INPUT_DIR2} \      
      ${OUTPUT_DIR} \
      ${CHR_LIST} \
      ${CATS}
  "


  # Cluster array execution
  JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}

#done
 
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
