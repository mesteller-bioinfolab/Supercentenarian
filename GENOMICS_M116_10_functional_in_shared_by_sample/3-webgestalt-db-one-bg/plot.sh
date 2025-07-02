

#DB=$(ls real/outputs | grep damaging | cut -d '_' -f2 | cut -d '.' -f1)
#NUM_DB=$(ls real/outputs | grep damaging | wc -l)

# per database  a facet


module load R

Rscript ./real/scripts/plot.R

