# sc17

#GENESET=papers_sc17
GENESET=all_longevity

CATS=( all coding damaging moderate modifier_cd modifier_nc altering cadd15 ncnc_cadd15 )

for CAT in ${CATS[@]}
do

cat real/outputs2/2TISSUES.M116.${GENESET}.overlaps.table | grep -w ${CAT} | cut -d ' ' -f5 > m116.${CAT}_${GENESET}

cat ../../../1000GIBS/10-functional/1-check-overlap/real/outputs2/1000GIBS.*.${GENESET}.overlaps.table | grep -w ${CAT} | cut -d ' ' -f5 > 1000g.${CAT}_${GENESET}

done







# all longevity (manel)

CATS=( damaging moderate modifier_cd modifier_nc altering cadd15 ncnc_cadd15 )

for CAT in ${CATS[@]}
do

cat real/outputs2/2TISSUES.M116.all_longevity.overlaps.table | grep -w ${CAT} | cut -d ' ' -f5 > m116_${CAT}_long

cat ../../../1000GIBS/10-functional/1-check-overlap/real/outputs2/1000GIBS.*.all_longevity.overlaps.table | grep -w ${CAT} | cut -d ' ' -f5 > 1000g.${CAT}_long

done

