#!/bin/bash
#
#SBATCH -p normal # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J 1000G
#SBATCH --mem 200G # memory pool for all cores
#SBATCH -t 9-23:59# time (D-HH:MM)
#SBATCH -o ./log.%j.out # STDOUT
#SBATCH -e ./log.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


# Get VCFs and indexes
FTP_SITE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot

for CHR in $(seq 1 22)
do
   FILE=$FTP_SITE/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants.annotated.vcf.gz
   wget $FILE 
   wget $FILE.tbi
   sleep 6
done


# X
wget $FTP_SITE/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.annotated.vcf.gz
wget $FTP_SITE/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.recalibrated_variants.annotated.vcf.gz.tbi




exit

# MT no hay

#-------
#module load PLINK/2.0
#cat ../individuals.txt | grep IBS | grep 30x | grep female | grep -v MSL | cut -f1 | sort -u > ../IBS_female.id 
 
#for CHR in $(seq 1 23 | sed 's/23/X/')
do
   FILE=$FTP_SITE/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.recalibrated_variants.annotated.vcf.gz
   plink2 --vcf $FILE --keep-fam IBS_female.id --geno 0.05 --make-bed --out tmp
   plink2 --bfile tmp --freq --out EUR_${CHR}.frq 
   #rm EUR_$CHR.{log,nosex} tmp.*
#done
 
awk 'FNR==1 && NR!=1{next;}{print}' *.frq | cut -f1-4 > EUR.freq
#rm EUR_*.frq





