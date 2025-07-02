rm -rf ./real/outputs3/*
 

# Run fisher test per individual, compara category voi with its recefence biotype
# Might as well filter this reference biotype by cat voi bg (seen in other 1000G in such VOI cat)


INPUT_DIR=./real/outputs2
OUTPUT_DIR=./real/outputs3



GENELISTS_DIR=./real/tmp

CAT_BG_DIR=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/outputs

cat ./real/inputs/cats.txt | grep -v "coding_with_nc" | grep -v "non_coding" | grep -v "all" > ./real/inputs/voi_cats.txt

CATS=./real/inputs/voi_cats.txt


TISSUE="1000GIBS"
ID_LIST=real/inputs/id.list

BACKGROUD_TYPE="combined"



#echo TISSUE IDD GENELIST CATEGORY PVALUE OR CI_L CI_U  > ${OUTPUT_DIR}/fet.table.header
echo TISSUE IDD GENELIST CATEGORY TOT_BG IN_BG TOT_CAT IN_CAT PVALUE OR CI_L CI_U  > ${OUTPUT_DIR}/fet.table.header


CATEGORIES=$(cat ${CATS})


for CAT in ${CATEGORIES[@]}

do



# COMMAND

COMMAND=" \
	./real/scripts/run_fet.sh \
      ${INPUT_DIR} \
      ${OUTPUT_DIR} \
      ${TISSUE} \
      ${ID_LIST} \
      ${GENELISTS_DIR} \
      ${CAT} \
      ${BACKGROUD_TYPE}
     "

	
  # Cluster array execution
  JOBS_COUNT=$(cat ${ID_LIST} | wc -l)
  sbatch --array=1-${JOBS_COUNT} ${COMMAND}
  





done
