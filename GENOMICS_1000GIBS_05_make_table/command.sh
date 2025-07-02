#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  rm -fR ./real/outputs/*


  TISSUE=1000GIBS

  # Configuration
  
  #module load BCFtools/1.9-foss-2018b  
  #bcftools query -f '%CHROM\n' ../0-data/${TISSUE}/M116.GATK.snp.vcf.gz  | sort -u > ./real/inputs/chrs.list
  
  module load GATK/4.1.6.0-GCCcore-8.3.0-Java-11


  #cat ../03-vep-annot/real/inputs/chrs.list | grep -v chrM | grep -v chrX | grep -v chrY > ./real/inputs/chrs.list

   cat ../03-vep-annot/real/inputs/chrs.list > ./real/inputs/chrs.list
   CHR_LIST=./real/inputs/chrs.list
   

  INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/04-cadd-annot/real/outputs
  OUTPUT_DIR=./real/outputs
  
  
   IDS=$(ls ${INPUT_DIR}/${TISSUE})

 
  # Execution 
  mkdir ${OUTPUT_DIR}/${TISSUE}


  # one job per ID in tissue
  
for ID in ${IDS}
do
echo $ID

mkdir ${OUTPUT_DIR}/${TISSUE}/${ID}


  # COMMAND
  COMMAND=" \
    ./real/scripts/vep2table.sh \
      ${CHR_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID}
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
