#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -rf ./real/inputs/*
  rm -rf ./real/tmp/*
  rm -fR ./real/outputs/*


  TISSUE=blood

  # Configuration
  
  #module load BCFtools/1.9-foss-2018b  
  #bcftools query -f '%CHROM\n' ../0-data/${TISSUE}/M116.GATK.snp.vcf.gz  | sort -u > ./real/inputs/chrs.list


  CHR_LIST=../02-vep-annot/real/inputs/chrs.list

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/08-filter-category-summary/real/outputs
  OUTPUT_DIR=./real/inputs
  
  
   IDS=M116
 
  # Execution 
  mkdir ${OUTPUT_DIR}




  # COMMAND
  COMMAND=" \
    ./real/scripts/count_samples.sh \
      ${CHR_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${IDS}
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
