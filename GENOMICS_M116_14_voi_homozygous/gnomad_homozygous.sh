#cat ../../gnomad_genomes_v4/female.homozygous.chr*  

#cat ../../gnomad_genomes_v4/female.homozygous.* | head | tr ' ' '_'  > real/inputs/gnomad.homozygous


# verify existence homozygous in gnomad v4

#cat real/outputs/2TISSUES.M116.ALTERING.homozygous.variants  | grep -f real/inputs/gnomad.homozygous 

# chr with homo

CHRS=( chr2 chr3 chr5 chr6 chr11 )
for CHR in ${CHRS[@]}
do
cat ../../gnomad_genomes_v4/female.homozygous.${CHR} | tr ' ' '_'  > real/inputs/gnomad.homozygous.${CHR}
#cat ../../gnomad_genomes_v4/female.homozygous.chr11 | tr ' ' '_'  > real/inputs/gnomad.homozygous.chr11
done
