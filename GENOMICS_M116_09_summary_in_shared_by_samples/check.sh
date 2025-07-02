# Genes -----

CATEGORY=altering

# common to 2 tissues
cat real/outputs/*/M116/*.rare_${CATEGORY}.genevar | cut -f1 | grep -v -w ENSG | sort -u | wc -l

# vs
TISSUED=blood
cat real/outputs/${TISSUED}/M116/*.rare_${CATEGORY}.genevar | cut -f1 | grep -v -w ENSG | sort -u | wc -l

TISSUED=saliva

TISSUED=urine
