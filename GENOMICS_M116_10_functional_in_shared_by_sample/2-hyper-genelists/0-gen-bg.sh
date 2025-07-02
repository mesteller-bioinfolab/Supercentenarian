
ls ../../../1000GIBS/08-filter-category-summary/real/outputs/1000GIBS/ > ./real/inputs/1000g.ids.txt



ID_LIST=./real/inputs/1000g.ids.txt


JOBS_COUNT=$(cat ${ID_LIST} | wc -l)

COMMAND="./real/scripts/gen.coding.bg_split.sh ${ID_LIST}"


eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
