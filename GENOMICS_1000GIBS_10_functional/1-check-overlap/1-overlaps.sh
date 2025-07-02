
#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  #rm -fR ./real/outputs/*
  
  # Genelists
 
  
  cp $( ls ../0-gene-lists/real/outputs/manel/*  ) ./real/tmp
  cp $( ls ../0-gene-lists/real/outputs/papers/* | grep -v chinese ) ./real/tmp
  cat ./real/tmp/*_ensg.txt | sort -u > ./real/tmp/all_longevity_ensg.txt  


  GENELISTS_DIR=./real/tmp


  TISSUE=1000GIBS



  CHR_LIST=../../02-vep-annot/real/inputs/chrs.list

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/09-summary/real/outputs


  # esto esta hecho ya en 10 x
  #AUTOSOMAL_GENES_FILE=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/summary_tables


  OUTPUT_DIR=./real/outputs

    
  #ID=M116

  ls -l ../../09-summary/real/outputs/ | grep txt | cut -d '.' -f2 | sort -u > ./real/inputs/id.list
  ID_LIST=./real/inputs/id.list

  CAT_FILE=./real/inputs/cats.txt

 
  # Execution 
  mkdir ${OUTPUT_DIR}


  # COMMAND
  COMMAND=" \
    ./real/scripts/overlap.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID_LIST} \
      ${GENELISTS_DIR} \
      ${CAT_FILE} \
      ${AUTOSOMAL_GENES_FILE}
  "

  # Cluster array execution
  JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  exit


  # Cluster execution
  eval sbatch ${COMMAND}
  exit             

} 


# FUNCTIONS ==========================================================

main
