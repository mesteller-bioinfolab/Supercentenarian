
#rm -rf ./real/outputs/*

INPUTS=/gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/12-rare-var-analyses/09-overlap-plot


#ls ../09-overlap-plot/*.inmain.genes | cut -d "_" -f2,3 | cut -d '.' -f1 > ./real/inputs/cats.txt
CAT_LIST=./real/inputs/cats.txt
SET_LIST=./real/inputs/sets.txt


COMMAND="./real/scripts/resampling_res_cat.sh ${INPUTS} ${CAT_LIST} ${SET_LIST}"


   #eg Rscript resampling_cast_cat.R inmain

# Cluster array execution
JOBS_COUNT=$(cat ${CAT_LIST} | wc -l)
eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
