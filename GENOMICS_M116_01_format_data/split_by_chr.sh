module load BCFtools/1.9-foss-2018b

TISSUE=$1
ID=$2
OUTDIR=./real/outputs

mkdir ${OUTDIR}/${TISSUE}
mkdir ${OUTDIR}/${TISSUE}/${ID}

VCF_IN=../00-data/wgs_data/${TISSUE}/${ID}.GATK.snp.vcf.gz

VCF_OUT_STEM=${OUTDIR}/${TISSUE}/${ID}/${TISSUE}.${ID}.by_chr

# or alll in chr.list


CHRS=$(bcftools query -f '%CHROM\n' ${VCF_IN} | sort -u)

#for i in {1..22}
for i in ${CHRS[@]}
do
#bcftools view ${VCF_IN} --regions chr${i} -o ${VCF_OUT_STEM}_chr${i}.vcf.gz -Oz
bcftools view ${VCF_IN} --regions ${i} -o ${VCF_OUT_STEM}_${i}.vcf.gz -Oz
done

