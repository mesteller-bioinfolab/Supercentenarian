
#2.9 years ago
wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz ./
zcat Homo_sapiens.GRCh38.104.chr.gtf.gz | grep "gene_biotype \"protein_coding\"" > Homo_sapiens.GRCh38.104.protein_coding.gtf


awk '$3=="gene"{print $10}' Homo_sapiens.GRCh38.104.protein_coding.gtf | sed 's/\"//g' | sed 's/\;//g'| sort -u  > Homo_sapiens.GRCh38.104.protein_coding.txt

