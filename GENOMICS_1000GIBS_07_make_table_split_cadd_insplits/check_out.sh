OUTS=$(ls real/outputs/1000GIBS | cut -d'.' -f2 | sort -u)

echo ID CHRS > check.txt
for O in ${OUTS}
do

echo ${O} $(ls real/outputs/1000GIBS/${O} | grep split2 | wc -l) >> check.txt

la real/outputs/1000GIBS/${O} | grep split2  >> ${O}.sizela.txt

done


25332341_4

25332336_17

25332310_3 

25332203_18
25332203_21
HG01626   5,8

25332183_19


ID CHRS
HG01519 21: 6
HG01626 20: 5,8
HG01685 21: 11
HG02222 21: 4
HG02232 21: 12
