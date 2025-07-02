module load VEP/102.0-foss-2019a-Perl-5.28.1 


# EUR AF filters



# "EUR_AF < 0.01 or not EUR_AF"
filter_vep -i ../1-vep-annot/real/outputs/blood/M116/blood.M116.chr1.vep -o filtered.blood.M116.g1000.chr1 --force_overwrite --filter "EUR_AF < 0.01 or not EUR_AF" --format vcf --vcf_info_field ANN --only_matched


# "gnomAD_NFE_AF < 0.01 or not gnomAD_NFE_AF"
filter_vep -i ../1-vep-annot/real/outputs/blood/M116/blood.M116.chr1.vep -o filtered.blood.M116.gnomad.chr1 --force_overwrite --filter "gnomAD_NFE_AF < 0.01 or not gnomAD_NFE_AF" --format vcf --vcf_info_field ANN --only_matched



# gnomAD_NFE_AF < 0.01 or EUR_AF < 0.01 or (not gnomAD_NFE_AF and not EUR_AF)
filter_vep -i ../1-vep-annot/real/outputs/blood/M116/blood.M116.chr1.vep -o filtered.blood.M116.eur.chr1 --force_overwrite --filter "EUR_AF < 0.01 or gnomAD_NFE_AF < 0.01 or (not gnomAD_NFE_AF and not EUR_AF)" --format vcf --vcf_info_field ANN --only_matched
