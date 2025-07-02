#!/bin/bash

# Main Code --------------------------------------------------------

main() {

  #rm -rf ./real/outputs/*


  # Configuration
    
  #module load GATK/4.1.6.0-GCCcore-8.3.0-Java-11

  CHR_LIST=./real/inputs/chrs.list
  
  # no, did not contain all it seems
  #INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/00-data/CCDG_14151_B01_GRM_WGS_2020-08-05_unnannotated
   INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/00-data/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_recalibrated_variants.annotated

  OUTPUT_DIR=./real/outputs

  TMP_DIR=./real/tmp


  
  # Execution 

#./real/scripts/filter_af_keepcommon.sh \

  # COMMAND
  COMMAND=" \
      ./real/scripts/common_snvs_annotated.sh \    
      ${CHR_LIST} \
      ${INPUT_DIR} \
      ${TMP_DIR} \
      ${OUTPUT_DIR}
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
