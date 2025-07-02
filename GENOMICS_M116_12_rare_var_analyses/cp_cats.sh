#ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 >  ./*/real/inputs/cats.txt

ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > 00-pergene-totals-vectors/real/inputs/cats.txt
ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > 01-CAST/real/inputs/cats.txt
ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > 02-MZ/real/inputs/cats.txt
ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > 03-rank-test/real/inputs/cats.txt
ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > 04-one-sample-prop/real/inputs/cats.txt


# also to other folders that need it:
# genelists
# papers
# webgestalt
ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/2-hyper-genelists/real/inputs/versions.txt
ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/2-hyper-papers/real/inputs/versions.txt
ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/3-webgestalt-db-gmt/real/inputs/versions.txt
ls /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/09-summary-in-shared-by-samples/real/outputs/blood/M116/blood.M116.chr12.* | grep dput | cut -d '.' -f4 | cut -d '_' -f3,4 > /gpfs42/projects/lab_anavarro/disease_pleiotropies/supercentenarian/Olot/10-functional-in-shared-by-sample/x-goof_red_genelists/real/inputs/cats.txt

