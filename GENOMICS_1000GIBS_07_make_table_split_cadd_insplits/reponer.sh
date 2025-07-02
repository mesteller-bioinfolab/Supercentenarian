IDS=$(cat pendingsplit | cut -f2 | sort -u )    #  grep -v HG01501)    # MANUALLY ADD ALL BUT CHR14 THATS RUNNING AGAIN NO REASON

IDS=$(cat pendingsplit | cut -f2 | sort -u | grep HG01501)

for IDD in ${IDS}
do

no | cp real/tmp/1000GIBS/${IDD}/1000GIBS.${IDD}.*.insplits real/tmp_1308_2/1000GIBS/${IDD}

done
