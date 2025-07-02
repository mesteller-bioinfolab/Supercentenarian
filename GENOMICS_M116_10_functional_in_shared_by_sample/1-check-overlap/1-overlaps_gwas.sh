
#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  #rm -fR ./real/tmp/*
  
  # Genelists
  #cp $( ls ../0-gene-lists/real/outputs/manel/*  ) ./real/tmp
  #cp $( ls ../0-gene-lists/real/outputs/papers/* | grep -v chinese ) ./real/tmp
  #cat ./real/tmp/*_ensg.txt | sort -u > ./real/tmp/all_longevity_ensg.txt  
  
  cp ../0-gene-lists/real/outputs/gwas/gwas_genes_ensg.txt ./real/tmp	


  GENELISTS_DIR=./real/tmp
  TISSUE="2TISSUES"



  CHR_LIST=../../02-vep-annot/real/inputs/chrs.list

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/summary_tables

  # esto esta hecho ya en 10 x
  #AUTOSOMAL_GENES_FILE=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/summary_tables


  OUTPUT_DIR=./real/outputs
    
  ID=M116
 
  CAT_FILE=./real/inputs/othercats.txt

 
  # Execution 
  mkdir ${OUTPUT_DIR}


  # COMMAND
  COMMAND=" \
    ./real/scripts/overlap2.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID} \
      ${GENELISTS_DIR} \
      ${CAT_FILE} \
      ${AUTOSOMAL_GENES_FILE}
  "

  # Cluster array execution
  #JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
  #eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  #exit


  # Cluster execution
  eval sbatch ${COMMAND}
  exit             

} 


# FUNCTIONS ==========================================================

main
