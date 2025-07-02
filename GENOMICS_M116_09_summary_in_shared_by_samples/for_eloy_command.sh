#mkdir real/outputs_fixedfulltables
#TMP_DIR=real/outputs_fixedfulltables

cat ../02-vep-annot/real/inputs/chrs.list | grep -v chrY > real/inputs/chrs.list


CHROMS=( $(cat real/inputs/chrs.list) )


for CHROM in ${CHROMS[@]}

do
#cat real/outputs/blood/M116/blood.M116.chr22.eur_rare_all.fulltable  | head -n1 > 2TISSUES.${CHROM}.fulltable



COMMAND="for_eloy.sh ${CHROM}"

sbatch ${COMMAND}

done
