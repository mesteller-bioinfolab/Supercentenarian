ID_LIST=( $(cat ./real/inputs/ids_list.txt) )


for IDD in ${ID_LIST[@]}
do
cat real/inputs/M116.chr20.genevar | head -n1 > ${IDD}.genevars
cat real/inputs/${IDD}.*.genevar | grep -v  -w symbol >> ${IDD}.genevars
done
