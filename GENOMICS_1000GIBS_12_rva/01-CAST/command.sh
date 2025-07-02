#!/bin/bash

# Main Code --------------------------------------------------------

main(){
  
  # Creates and cleans directory structure
  clean_directory
  
  #cp /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/11-g-matrix/real/outputs/rare.gmat /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/11-g-matrix/real/outputs/all.gmat


  # Configuration

  TISSUE=1000GIBS
  ls -l ../../09-summary/real/outputs/ | grep txt | cut -d '.' -f2 | sort -u > ./real/inputs/id.list
  ID_LIST=./real/inputs/id.list


  CAT_FILE=./real/inputs/cats.txt
  INPUTS_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/11-g-matrix/real/outputs
  OUTPUTS_DIR=./real/outputs


  # COMMAND
  COMMAND=" \
    ./real/scripts/run_cast.sh \
    ${INPUTS_DIR} \
    ${OUTPUTS_DIR} \
    ${TISSUE} \
    ${ID_LIST} \
    ${CAT_FILE}
  "
  
   
  # Execution 
  # Cluster array execution
  JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  #exit

  # Direct execution
  #eval bash ${COMMAND}
  #exit 

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
