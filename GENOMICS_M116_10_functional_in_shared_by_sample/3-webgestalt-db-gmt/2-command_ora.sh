
#!/bin/bash

# Main Code --------------------------------------------------------

main(){

  rm -rf real/tmp/*
  rm -rf real/outputs/*

  # Creates and cleans directory structure
  #clean_directory



  # Configuration

   VERSIONS_LIST=./real/inputs/versions.txt

  
   
   INPUTS_DIR=../../09-summary-in-shared-by-samples/real/summary_tables


   # previously generated here
   GENESETS_DIR=./real/inputs
   OUTPUTS_DIR=./real/outputs


  
   
  #CAT_BG_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/outputs
  CAT_BG_DIR=../1-check-overlap/real/outputs/compound_bg

  # COMMAND
  COMMAND=" \
    ./real/scripts/run_ora.sh \
    ${INPUTS_DIR} \
    ${GENESETS_DIR} \
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
