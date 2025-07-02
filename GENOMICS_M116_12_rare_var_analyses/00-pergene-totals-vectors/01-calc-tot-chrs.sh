#!/bin/bash

# Main Code --------------------------------------------------------

main(){
  
  # Creates and cleans directory structure
  clean_directory

  ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > ./real/inputs/cats.txt


  # Configuration
  CATEGORIES=$( cat ./real/inputs/cats.txt )

  #CATEGORIES=( damaging moderate modifier_cd modifier_nc altering low low_only )
  INPUTS_DIR=./real/inputs
  OUTPUTS_DIR=./real/tmp
  
  cat ../../02-vep-annot/real/inputs/chrs.list | grep -v chrM | grep -v chrY | grep -v chrX > ./real/inputs/chrs.list
  CHR_LIST=./real/inputs/chrs.list

  for CATEGORY in ${CATEGORIES[@]}
  do
  
  # COMMAND
  COMMAND=" \
    ./real/scripts/tot.sh \
    ${INPUTS_DIR} \
    ${OUTPUTS_DIR} \
    ${CATEGORY} \
    ${CHR_LIST}
  "
  
  # Execution 
  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  
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
