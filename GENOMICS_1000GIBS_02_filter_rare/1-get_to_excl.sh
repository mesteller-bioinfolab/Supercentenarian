# get AF_EUR_unrel > 0.015 list in the format of these VCF ID columsn to filter
# with bcftools view -e 'ID=@excl.snp' real/outputs/1000GIBS/HG02222/1000GIBS.HG02222.chr1.filtered.gz
# this instead of using own VCF column AF_EUR_unrel for no reason, homogeneity with M116 treatment, 
# in case the AF varies among versions , here or annotated files in
# ../../1000GIBS/00-data/common_snvs/real/outputs
#cat  ../../1000GIBS/00-data/common_snvs/real/outputs/1000G.chr6.eur.common.table | head | cut -d ' ' -f1,2,3,4 |  tr ' ' ':' | awk '{ gsub("chr", "") ; print $0 }'
 


INPUTS_DIR=../../1000GIBS/00-data/common_snvs/real/outputs
OUTPUTS_DIR=./real/inputs
CHR_LIST=./real/inputs/chrs.list

COMMAND="./real/scripts/get_to_excl.sh ${INPUTS_DIR} ${OUTPUTS_DIR} ${CHR_LIST}"

JOBS_COUNT=$(cat ${CHR_LIST} | wc -l)
eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}

