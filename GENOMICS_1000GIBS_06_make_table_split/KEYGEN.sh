# ANN field that has 28 fields (KEY2)

To use python parser for VEP ANN field in gatk VariantTo Table formatted file:
- Manually copied field names in VCF ANN field to a file named INFO_ANN with nano or:
cat ../5-cadd-annot/real/outputs/1000GIBS/HG01504/1000GIBS.HG01504.chr12.vep | grep INFO | grep ANN | cut -d ':' -f2 | cut -d ' ' -f2 | cut -d '"' -f1 > ./real/inputs/INFO_ANN
- Changed delimiter: cat ./real/inputs/INFO_ANN | tr '|' '\n'  > ./real/inputs/CSQ_KEY_
- Changed CSQ by ANN in the script
- Run script using CSQ_KEY_


