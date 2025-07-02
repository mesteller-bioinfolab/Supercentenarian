#!/bin/bash

# Main Code --------------------------------------------------------

main() {
  #rm -fR ./real/outputs/*
 
  ANN_KEY=./real/inputs/CSQ_KEY_

  TISSUE=1000GIBS

  # Configuration
  
  #module load BCFtools/1.9-foss-2018b  
  #bcftools query -f '%CHROM\n' ../0-data/${TISSUE}/M116.GATK.snp.vcf.gz  | sort -u > ./real/inputs/chrs.list
  
  #module load GATK/4.1.6.0-GCCcore-8.3.0-Java-11

  # CHR subset
  #cat ../4-vep-annot/real/inputs/chrs.list | grep chr21 > real/inputs/chrs.list
  #cat ../4-vep-annot/real/inputs/chrs.list | grep chr22 >> real/inputs/chrs.list
  #CHR_LIST=real/inputs/chrs.list

  # All CHRs
  #cat ../4-vep-annot/real/inputs/chrs.list | grep -v chrM | grep -v chrX | grep -v chrY > ./real/inputs/chrs.list2
  CHR_LIST=./real/inputs/chrs.list

  #INPUT_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/1000GIBS/7-make-table-split/real/outputs
   INPUT_DIR=./real/tmp
   OUTPUT_DIR=./real/tmp
  
  
   #IDS=HG01510
   IDS=$(ls ${INPUT_DIR}/${TISSUE} | grep -v HG02230 | grep -v HG01507 )
   #IDS=HG02230
   
  # Execution 
  #mkdir ${OUTPUT_DIR}/${TISSUE}
  #rm -fR ./real/outputs/${TISSUE}


  # one job per ID in tissue
  
for ID in ${IDS}
do
echo $ID

#mkdir ${OUTPUT_DIR}/${TISSUE}/${ID}

CHR_SPLIT_LIST=./real/inputs/${TISSUE}.${ID}.chrs_split.list

  # COMMAND
  COMMAND=" \
    ./real/scripts/ann-split-insplits-lite.sh \
      ${CHR_SPLIT_LIST} \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID} \
      ${ANN_KEY}
  "

  # Cluster array execution
  JOBS_COUNT=$(cat ${CHR_SPLIT_LIST} | wc -l)
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
