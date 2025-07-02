cat real/inputs/cats.txt | grep -v low | grep -v coding > ./real/inputs/voi_cats.txt



CATEGORIES=$(cat ./real/inputs/voi_cats.txt)


echo TOTAL M116 > ${CATEGORY}_voi_tot_m116.txt

for CATEGORY in ${CATEGORIES[@]}


do


CAT_DIFF_GENES=( $(cat /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/12-rare-var-analyses/09-overlap-plot/real/outputs/a5_${CATEGORY}.inmain.genes) )


for GENE in ${CAT_DIFF_GENES[@]}

do


echo ${GENE} $(cat /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/12-rare-var-analyses/00-pergene-totals-vectors/real/tmp/${CATEGORY}.*.tot_gene_vars | grep ${GENE} | cut -f2) $(cat /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/11-g-matrix/real/tmp/M116.*.id_${CATEGORY}_g  |  grep ${GENE} | cut -f2) >> ${CATEGORY}_voi_tot_m116.txt

done



done
