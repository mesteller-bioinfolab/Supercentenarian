#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -fR ./real/outputs/urine/*

  SUFF="fullsplit_AF2"

  TISSUE=urine

  # Configuration
  
  #module load BCFtools/1.9-foss-2018b  
  #bcftools query -f '%CHROM\n' ../0-data/${TISSUE}/M116.GATK.snp.vcf.gz  | sort -u > ./real/inputs/chrs.list


  CHR_LIST=../02-vep-annot/real/inputs/chrs.list

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/08-filter-category-summary/real/outputs
  TMP_DIR=./real/tmp
  OUTPUT_DIR=./real/outputs
  FILTER_DIR=../../1000GIBS/00-data/common_snvs/real/outputs  
  
   #IDS=M116
    IDS=($(ls ${INPUT_DIR}/${TISSUE}))

  # Execution 
  mkdir ${OUTPUT_DIR}
  mkdir ${OUTPUT_DIR}/${TISSUE}



  # one job per ID in tissue
  
for ID in ${IDS[@]}
do
echo ${ID}

mkdir ${OUTPUT_DIR}/${TISSUE}/${ID}

  # COMMAND
  COMMAND=" \
      ./real/scripts/crop-summary.sh \
      ${CHR_LIST} \
      ${INPUT_DIR} \
      ${TMP_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID} \
      ${SUFF} \
      ${FILTER_DIR}
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
