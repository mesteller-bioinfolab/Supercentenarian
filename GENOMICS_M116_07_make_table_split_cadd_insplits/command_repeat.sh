#IDD=M116
#TISSUED=saliva
#chr21_ad
#chr2_ap

IDD=T90
TISSUED=saliva
#chr21_ad
#chr2_aq
CHR=chr2_aq
KEY=./real/inputs/CSQ_KEY_

sbatch repeat_failed.sh ${IDD} ${TISSUED} ${CHR} ${KEY}
