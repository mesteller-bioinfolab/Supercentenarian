rm -rf ./real/outputs4/*


GENELISTS_DIR=./real/tmp

CATS=./real/inputs/cats.txt


OUTPUTS_DIR=./real/outputs4
mkdir ${OUTPUTS_DIR}





CATEGORIES=($(cat ${CATS}))

IBS_DIR=real/outputs2
M116_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/1-check-overlap_gwas/real/outputs2


#BORRARfor CATEGORY in ${CATEGORIES[@]}
#do
#ls ${IBS_DIR}/1000GIBS.*.${CATEGORY}.overlaps.table > ./real/inputs/1000GIBS.${CATEGORY}.overlaps.table
#done




GENELISTS=$(ls ${GENELISTS_DIR}/*.txt | cut -d '/' -f4 | rev | cut -d '_' -f2,3,4 | rev)


 
# Execution 
 

#ls ${IBS_DIR}/1000GIBS.*.${GENELIST}.overlaps.table > ./real/inputs/1000GIBS.${GENELIST}.overlaps.table





for GENELIST in ${GENELISTS[@]}
do

ls ${IBS_DIR}/1000GIBS.*.${GENELIST}.overlaps.table > ./real/inputs/1000GIBS.${GENELIST}.overlaps.tables.paths


  # COMMAND
  COMMAND=" \
    ./real/scripts/confidence.interval.sh \
      ${M116_DIR} \
      ${GENELIST} \
      ${CATS} \
      ${OUTPUTS_DIR}
  "

  sbatch ${COMMAND}
done
