#!/bin/bash

# Main Code --------------------------------------------------------

main() {


  # Configuration

  # CHRs 
  cat ../02-vep-annot/real/inputs/chrs.list | grep -v chrM | grep -v chrX | grep -v chrY > ./real/inputs/chrs.list
  

    CHR_LIST=./real/inputs/chrs.list 
   
  
   
   INPUT_DIR=./real/inputs
   OUTPUT_DIR=./real/outputs

 
    #echo M116 > ./real/inputs/ids_list.txt
    #ls -l ../../1000GIBS/09-summary/real/outputs/ | grep txt | grep -v -w 116 | cut -d ' ' -f14 | cut -d '.' -f2 >> ./real/inputs/ids_list.txt
   
    #cat ./real/inputs/ids_list.txt  | sed -n 2p  > ./real/inputs/ids_list1.txt
      
    #IDS_LIST=./real/inputs/ids_list.txt
  

  # Execution 
   
  # one job per CAT
  
#CATEGORIES=( damaging moderate modifier_cd modifier_nc altering low coding all )


CATEGORIES=$( cat ./real/inputs/cats.txt ) 

for CATEGORY in ${CATEGORIES[@]}
do
echo ${CATEGORY}

  # COMMAND
  COMMAND=" \
    ./real/scripts/join.chrs.cat.sh \
      ${CATEGORY} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR}
  "


  # Cluster array execution
  #JOBS_COUNT=$(cat ${IDS_LIST} | wc -l)
  #eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  sbatch ${COMMAND}

done
 
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
