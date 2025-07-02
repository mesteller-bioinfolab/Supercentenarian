TISSUE=saliva
ID=T90

# inputs
cat real/inputs/${TISSUE}.${ID}.chrs_split.list | cut -d'.' -f4 | sort -u > ${TISSUE}.${ID}.in
# done
la real/tmp/${TISSUE}/${ID}/${TISSUE}.${ID}.* | grep split2 | sort -u | cut -d'.' -f3 | sort -u >  ${TISSUE}.${ID}.done

# diff
comm -23 ${TISSUE}.${ID}.in ${TISSUE}.${ID}.done  >  ${TISSUE}.${ID}.todo

