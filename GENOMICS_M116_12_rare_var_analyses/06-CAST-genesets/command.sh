#!/bin/bash

# Main Code --------------------------------------------------------

main(){
  
  # Creates and cleans directory structure
  clean_directory
  
  # Configuration
  cp ../../10-functional-in-shared-by-sample/0-gene-lists/real/outputs/manel/*.txt  ./real/tmp 
  cp ../../10-functional-in-shared-by-sample/0-gene-lists/real/outputs/papers/*.txt  ./real/tmp
  cat ./real/inputs/*.txt  > ./real/tmp/longevity.txt
  cp ../../10-functional-in-shared-by-sample/0-gene-lists/real/outputs/gwas*.txt ./real/tmp


  VERSIONS_LIST=./real/inputs/cats.txt
  INPUTS_DIR=../01-CAST/real/outputs
  OUTPUTS_DIR=./real/outputs


  # COMMAND
  COMMAND=" \
    ./real/scripts/run_cast.sh \
    ${INPUTS_DIR} \
    ${OUTPUTS_DIR} \
    ${VERSIONS_LIST}
  "
  
   
  # Execution 
  # Cluster array execution
  JOBS_COUNT=$(cat ${VERSIONS_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  #exit

  # Direct execution
  eval bash ${COMMAND}
  exit 

  # Cluster execution
  #eval sbatch ${COMMAND}
  #exit             

} 



# FUNCTIONS ==========================================================


clean_directory(){
  # Creates and cleans directory structure
  
  if [ ! -d "real" ]; then

    mkdir real
    mkdir real/scripts
    mkdir real/inputs
    mkdir real/tmp
    mkdir real/outputs

  else

    if [ ! -d "real/scripts" ]; then
  mkdir real/scrips
    fi

    if [ ! -d "real/inputs" ]; then
  mkdir real/inputs
    fi

    if [ ! -d "real/tmp" ]; then
  mkdir real/tmp
    else
  rm -fR real/tmp/*
    fi

    if [ ! -d "real/outputs" ]; then
  mkdir real/outputs
    else
  rm -fR real/outputs/*
    fi

  fi
}
 
 
main
