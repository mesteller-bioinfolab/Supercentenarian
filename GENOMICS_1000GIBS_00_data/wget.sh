#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J 1000G
#SBATCH --mem 200G # memory pool for all cores
#SBATCH -t 9-23:59# time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address






# Get VCFs

#FTP_SITE=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502

#FTP_SITE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
FTP_SITE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased

#wget $FTP_SITE/integrated_call_samples_v3.20130502.ALL.panel


for CHR in $(seq 1 22)
do
   #FILE=$FTP_SITE/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
   FILE=$FTP_SITE/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
   wget $FILE $FILE.tbi
   sleep 60
done


# X
wget $FTP_SITE/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.shapeit2-duohmm-phased.vcf.gz


#-------
module load PLINK/2.0

# Identify EUR and calculate AF. (ignoring variants with call rate < 95% and NO: that fail HW E at p<1e-06). PLINK1.9 or 2
# no --hwe 1e-6
grep EUR integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > EUR.id  ## 503 Europeans
 
for CHR in $(seq 1 23 | sed 's/23/X/')
do
   #FILE=ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
   FILE=CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
   plink2 --vcf $FILE --keep-fam EUR.id --geno 0.05 --make-bed --out tmp
   plink2 --bfile tmp --freq --out EUR_${CHR}.frq 
   rm EUR_$CHR.{log,nosex} tmp.*
done
 
awk 'FNR==1 && NR!=1{next;}{print}' *.frq | cut -f1-4 > EUR.freq
#rm EUR_*.frq





