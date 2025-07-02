

for

cat header  > ./real/inputs/2TISSUES.M116.${CHR}.eur_rare_${CAT}.fulltable
cat ../09-summary-in-shared-by-samples/real/outputs/*/M116/*.${CHR}.eur_rare_${CAT}.fulltable | sort -u   >>  ./real/inputs/2TISSUES.M116.${CHR}.eur_rare_${CAT}.fulltable
