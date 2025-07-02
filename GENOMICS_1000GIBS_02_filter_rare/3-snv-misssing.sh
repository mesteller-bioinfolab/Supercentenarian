# get AF_EUR_unrel > 0.015 list in the format of these VCF ID columsn to filter
# with bcftools view -e 'ID=@excl.snp' real/outputs/1000GIBS/HG02222/1000GIBS.HG02222.chr1.filtered.gz
# this instead of using own VCF column AF_EUR_unrel for no reason, homogeneity with M116 treatment, 
# in case the AF varies among versions , here or annotated files in
# ../../1000GIBS/00-data/common_snvs/real/outputs
#cat  ../../1000GIBS/00-data/common_snvs/real/outputs/1000G.chr6.eur.common.table | head | cut -d ' ' -f1,2,3,4 |  tr ' ' ':' | awk '{ gsub("chr", "") ; print $0 }'
 


INPUTS_DIR=../01-format-data/real/outputs
TMPS_DIR=./real/inputs
OUTPUTS_DIR=./real/outputs
TISSUE=1000GIBS

CHR_ID_LIST=pending


#IDS=$(ls ${INPUTS_DIR}/${TISSUE})
# Im unable to specify per chr snp file in 
#so have to load the full list every time
#yes | rm ./real/inputs/1000G.common.snps
#cat real/inputs/*.common.id > ./real/inputs/1000G.common.snps
#cat real/inputs/*.common.id | grep -v "CHROM:POS:REF:ALT" > snpsexcl
# Execution 
#mkdir ${OUTPUTS_DIR}
#mkdir ${OUTPUTS_DIR}/${TISSUE}
# one job per ID in tissue
#for ID in ${IDS}
#do
#echo $ID
#mkdir ${OUTPUTS_DIR}/${TISSUE}/${ID}



COMMAND="./real/scripts/snvs_missing.sh ${INPUTS_DIR} ${TMPS_DIR} ${OUTPUTS_DIR} ${CHR_ID_LIST} ${TISSUE}"


JOBS_COUNT=$(cat ${CHR_ID_LIST} | wc -l)
eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}
