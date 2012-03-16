lftp -c open -e "mirror -e -i ffn -i faa -i gbk -i ptt -x scaffold -n /genomes/Bacteria_DRAFT/ /work/rec3141/gene_clustering/all_genomes/contig" ftp://ftp.ncbi.nih.gov

lftp -c open -n /genomes/Bacteria_DRAFT/ -e "find ./ | grep ffn.tgz | cut -f3 -d'/' > allcontigs.txt" ftp://ftp.ncbi.nih.gov
