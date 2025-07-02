
#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  
  rm -fR ./real/outputs2/*
  

  # Genelist alreday there
  #cp $( ls ../0-gene-lists/real/outputs/*/* | grep -v chinese ) ./real/tmp
  #cat ../0-gene-lists/real/outputs/*/* | sort -u > ./real/tmp/all_longevity_ensg.txt  



  GENELISTS_DIR=./real/tmp
  TISSUE=1000GIBS



  #CHR_LIST=../../02-vep-annot/real/inputs/chrs.list



  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/09-summary/real/outputs


  # esto esta hecho ya en 10 x
  #AUTOSOMAL_GENES_FILE=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/summary_tables


  OVER_DIR=./real/outputs
    
 
  ls -l ../../09-summary/real/outputs/ | grep txt | cut -d '.' -f2 | sort -u > ./real/inputs/id.list
  ID_LIST=./real/inputs/id.list

  CAT_FILE=./real/inputs/cats.txt


  OUTPUT_DIR=./real/outputs2
  mkdir ${OUTPUT_DIR}


  echo GENESET GENESET_SIZE CATEGORY TOT_CAT IN_CAT OUT_CAT  > ${OUTPUT_DIR}/overlaps.table.header


  #echo GENESET CATEGORY GENESET_SIZE CAT_TOT CAT_IN CAT_OUT  > overlaps.table



  GENELISTS=$(ls ${GENELISTS_DIR}/*.txt | cut -d '/' -f4 | rev | cut -d '_' -f2,3,4 | rev)


 
  # Execution 
 

for GENELIST in ${GENELISTS[@]}
do

  # COMMAND
  COMMAND=" \
    ./real/scripts/explorer.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID_LIST} \
      ${GENELISTS_DIR} \
      ${GENELIST} \
      ${CAT_FILE} \
      ${OVER_DIR}
  "

  # Cluster array execution
  JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
  eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
 

done

exit

  # Cluster execution
  eval sbatch ${COMMAND}
  exit             

} 


# FUNCTIONS ==========================================================

main
