

CATS=$(cat ./real/inputs/cats.txt)
GE

for CAT in ${CATS}
do
comm -3 real/outputs/2TISSUES.M116.${CAT}.all_longevity.25chrs real/outputs/2TISSUES.M116.${CAT}.all_longevity > real/outputs/2TISSUES.M116.${CAT}.all.longevity.chrX

#comm -3 real/outputs/2TISSUES.M116.${CAT}.all_longevity.25chrs real/outputs/2TISSUES.M116.${CAT}.all_longevity > longevity.chrX
done
