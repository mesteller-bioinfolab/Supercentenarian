#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -fR ./real/outputs/*
 
  ANN_KEY=./real/inputs/CSQ_KEY_

  TISSUE=1000GIBS

  # Configuration
  
  #module load BCFtools/1.9-foss-2018b  
  #bcftools query -f '%CHROM\n' ../0-data/${TISSUE}/M116.GATK.snp.vcf.gz  | sort -u > ./real/inputs/chrs.list
  
  module load GATK/4.1.6.0-GCCcore-8.3.0-Java-11

  CHR_LIST=../03-vep-annot/real/inputs/chrs.list

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/05-make-table/real/outputs
  OUTPUT_DIR=./real/outputs
  
  
 
   IDS=$(ls ${INPUT_DIR}/${TISSUE})

 
  # Execution 
  mkdir ${OUTPUT_DIR}/${TISSUE}
  #rm -fR ./real/outputs/${TISSUE}


  # one job per ID in tissue
  
for ID in ${IDS}
do
echo $ID

mkdir ${OUTPUT_DIR}/${TISSUE}/${ID}


  # COMMAND
  COMMAND=" \
    ./real/scripts/ann-split.sh \
      ${CHR_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID} \
      ${ANN_KEY}
  "

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

main
