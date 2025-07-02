ls real/outputs/1000GIBS/ > ./real/inputs/ids.txt


IDS=$(cat ./real/inputs/ids.txt)

echo ID total > id_tables.txt



for ID in ${IDS}
do

echo $ID $(la real/outputs/1000GIBS/HG01501/ | grep table | wc -l) >> id_tables.txt

done
