
CHR_LIST=./real/inputs/chrs.list
INPUTS=./real/outputs
TMPS=./real/tmp2
OUTPUTS=./real/outputs2

TISSUE=1000GIBS

mkdir ${TMPS}
mkdir ${TMPS}/${TISSUE}


mkdir ${OUTPUTS}
mkdir ${OUTPUTS}/${TISSUE}




COMMAND="./real/scripts/get_to_exclude.sh ${CHR_LIST} ${INPUTS} ${TMPS} ${OUTPUTS}"


JOBS_COUNT=$(cat ${CHR_LIST} | wc -l)
eval sbatch --array=1-${JOBS_COUNT} ${COMMAND}




# exclude
#bcftools view  -i 'QD < 2.0 || FS > 60.0 || MQ < 40.0  || MQRankSum < -12.5 || ReadPosRankSum < -8.0' ${TISSUED}/${IDD}/${IDD}.${CHR}.vcf.gz > ${TISSUED}/${IDD}/${CHR}.to_exclude.vcf

# filter
#bcftools view  -e 'QD < 2.0 || FS > 60.0 || MQ < 40.0  || MQRankSum < -12.5 || ReadPosRankSum < -8.0' ${TISSUED}/${IDD}/${IDD}.${CHR}.vcf.gz  > ${TISSUED}/${IDD}/${IDD}.${CHR}.filtered.vcf.gz

# exclude table
#bcftools query -f '%CHROM %POS %REF %ALT\n' ${TISSUED}/${IDD}/${CHR}.to_exclude.vcf > ${TISSUED}/${IDD}/${CHR}.to_exclude.table

