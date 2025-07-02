# FILTER PASS
cat real/outputs/table.M116.chr1.good | head -n3 | awk '{ if ($7 == "PASS") print $0 }' 

# QUAL > 200
#cat real/outputs/table.M116.chr1.good | head -n3 | awk '{ if ($6 > 200) print $0 }' 



# ENSBL Gene ID 14
cat real/outputs/table.M116.chr1.good | head -n3 | cut -f14

# HGNC_ID 33
cat real/outputs/table.M116.chr1.good | head -n3 | cut -f33



# Feature_type 15

# Feature 16 

# Existing variation 27 
cat real/outputs/table.M116.chr1.good | head -n3 | cut -f27


# CANONICAL YES/""
cat real/outputs/table.M116.chr1.good | head -n3 | cut -f34

#ENSP 39

# SWISSPROT 40

# GENE_PHENO 44   either 1 or ""

