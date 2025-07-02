#!/bin/bash

# Main Code --------------------------------------------------------

main(){
  
  # Creates and cleans directory structure
  clean_directory


  # Configuration

  VERSIONS_LIST=./real/inputs/versions.txt
  INPUTS_DIR=../../09-summary-in-shared-by-samples/real/outputs_ya   
  OUTPUTS_DIR=./real/outputs


  # All coding genes in eur rare dataset
  #cat  ../../09-summary-in-shared-by-samples/real/outputs/*/M116/*.eur_rare_all.fulltable | grep Transcript | grep protein_coding | grep -v -w ID | cut -f16 | sort -u | grep -v -w "" 
  # or, same, precalculated
  #cat  ../../09-summary-in-shared-by-samples/real/outputs/*/M116/*.eur_rare_coding.genes | sort -u > ./real/inputs/M116_coding.genes.txt


  # All coding genes only in autosomes
  #cat  ../../09-summary-in-shared-by-samples/real/outputs/*/M116/*.eur_rare_all.fulltable | grep -v -w chrX | grep -v -w chrY | grep -v -w chrM | grep Transcript | grep protein_coding | grep -v -w ID | cut -f16 | sort -u | grep -v -w ""  >  ./real/inputs/M116_coding.genes.txt
  #BACKGROUND=./real/inputs/M116_coding.genes.txt


  CAT_BG_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/outputs


  # COMMAND
  COMMAND=" \
    ./real/scripts/run.sh \
    ${INPUTS_DIR} \
    ${OUTPUTS_DIR} \
    ${VERSIONS_LIST} \
    ${CAT_BG_DIR}
  "
  
   
  # Execution 
  # Cluster array execution
  JOBS_COUNT=$(cat ${VERSIONS_LIST} | wc -l)
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
