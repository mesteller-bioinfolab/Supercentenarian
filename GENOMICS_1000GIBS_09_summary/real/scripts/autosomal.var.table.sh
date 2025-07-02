INPUTS=./real/summary_tables


cat ${INPUTS}/eur_rare.variants | grep -v chrX | grep -v chrM > ${INPUTS}/eur_rare.autosomal.variants


cat ${INPUTS}/eur_rare_damaging.variants | grep -v chrX | grep -v chrM >  ${INPUTS}/eur_rare_damaging.autosomal.variants
cat ${INPUTS}/eur_rare_moderate.variants | grep -v chrX | grep -v chrM > ${INPUTS}/eur_rare_moderate.autosomal.variants
cat ${INPUTS}/eur_rare_modifier_cd.variants | grep -v chrX | grep -v chrM > ${INPUTS}/eur_rare_modifier_cd.autosomal.variants
cat ${INPUTS}/eur_rare_modifier_nc.variants | grep -v chrX | grep -v chrM > ${INPUTS}/eur_rare_modifier_nc.autosomal.variants
cat ${INPUTS}/eur_rare_altering.variants | grep -v chrX | grep -v chrM > ${INPUTS}/eur_rare_altering.autosomal.variants


