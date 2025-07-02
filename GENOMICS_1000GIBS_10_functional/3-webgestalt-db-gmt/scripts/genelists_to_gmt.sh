#!/bin/bash
#
#SBATCH -p haswell # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -J ora_genelists
#SBATCH --mem 16G # memory pool for all cores
#SBATCH -t 0-03:00 # time (D-HH:MM)
#SBATCH -o ./real/outputs/log.%j.out # STDOUT
#SBATCH -e ./real/outputs/log.%j.err # STDERR
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=claudia.vasallo@upf.edu # send-to address


#rm -rf ./real/inputs/*
#rm -rf ./real/tmp/*
#rm -rf ./real/outputs/*


module load R


# Databasaes in gmt 
#cp ../0-gene-lists/real/inputs/*.gmt real/inputs/
cp /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/13-functional-in-rva/3-webgestalt-db-gmt/real/inputs/*.gmt ./real/inputs


exit


# Genesets in files, to convert to gmt

GENELIST_DIR_MANEL=$(ls -d ../0-gene-lists/real/outputs/manel/*)

GENELIST_DIR_PAPERS=$(ls -d ../0-gene-lists/real/outputs/papers/*)



for G in ${GENELIST_DIR_MANEL[@]}
do
GN=$(basename $G)
echo $G
echo $GN
Rscript ./real/scripts/genelists_to_gmt.R ${G}


done 


for G in ${GENELIST_DIR_PAPERS[@]}
do
GN=$(basename $G)
echo $G
echo $GN
Rscript ./real/scripts/genelists_to_gmt.R ${G}
cat ./real/inputs/${GN}.gmt > ./real/inputs/${GN}.gmt

#cat msigdb.v2023.2.Hs.entrez.gmt >> ./real/inputs/${GN}.gmt
#cat ./real/inputs/c5.hpo.v2023.2.Hs.entrez.gmt >> ./real/inputs/${GN}.gmt

done


# BG genes
#G=real/inputs/BG_TWO_M116.genes
#GN=ref
#Rscript ./real/scripts/ref_to_gmt.R ${G} > ./real/inputs/${GN}.gmt



exit 

mkdir test


# ref plus longevity gmt
cat real/inputs/ref.gmt > real/inputs/ref_and_longevity.gmt
mv real/inputs/ref.gmt test
cat real/inputs/excel.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/excel.gmt test
cat real/inputs/hagr_genage_aging.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/hagr_genage_aging.gmt test
cat real/inputs/hagr_genage_human.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/hagr_genage_human.gmt test
cat real/inputs/hagr_cellage.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/hagr_cellage.gmt test
cat real/inputs/hagr_cellsignatures.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/hagr_cellsignatures.gmt test
cat real/inputs/longevitymap.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/longevitymap.gmt test
cat real/inputs/ngdc.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/ngdc.gmt test
cat real/inputs/papers_sc17.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/papers_sc17.gmt test
cat real/inputs/papers_chinese.gmt >> real/inputs/ref_and_longevity.gmt
mv real/inputs/papers_chinese.gmt test



