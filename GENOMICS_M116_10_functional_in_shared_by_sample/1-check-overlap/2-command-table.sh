
#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -fR ./real/outputs2/*
  
  # Genelists
 
  #cp $( ls ../0-gene-lists/real/outputs/*/* | grep -v chinese ) ./real/tmp
  #cat ../0-gene-lists/real/outputs/*/* | sort -u > ./real/tmp/all_longevity_ensg.txt  


  GENELISTS_DIR=./real/tmp
  TISSUE="2TISSUES"

  #CHR_LIST=../../02-vep-annot/real/inputs/chrs.list


  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/summary_tables


  # esto esta hecho ya en 10 x
  #AUTOSOMAL_GENES_FILE=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/summary_tables


  OVER_DIR=./real/outputs
    
 
  #ls -l ../../09-summary/real/outputs/ | grep txt | cut -d '.' -f2 | sort -u > ./real/inputs/id.list
  #ID_LIST=./real/inputs/id.list

   ID=M116
   CAT_FILE=./real/inputs/cats.txt


  OUTPUT_DIR=./real/outputs2
  mkdir ${OUTPUT_DIR}

  echo GENESET GENESET_SIZE CATEGORY TOT_CAT IN_CAT OUT_CAT  > ${OUTPUT_DIR}/overlaps.table.header

 
  # Execution 

   #GENELISTS=$(ls ${GENELISTS_DIR}/*.txt | cut -d '/' -f6 | cut -d '.' -f1)
    GENELISTS=$(ls ${GENELISTS_DIR}/*.txt | cut -d '/' -f4 | rev | cut -d '_' -f2,3,4 | rev) 


for GENELIST in ${GENELISTS[@]}
do

  # COMMAND
  COMMAND=" \
    ./real/scripts/explorer.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID} \
      ${GENELISTS_DIR} \
      ${GENELIST} \
      ${CAT_FILE} \
      ${OVER_DIR}
  "

  # Cluster array execution
  #JOBS_COUNT=$(cat ${CAT_FILE} | wc -l)
  #eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
 

  eval sbatch ${COMMAND}

done

exit


  # Cluster execution
  eval sbatch ${COMMAND}


  exit             

} 


# FUNCTIONS ==========================================================

main
