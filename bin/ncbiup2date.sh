#!/bin/bash

alldir=/work/rec3141/gene_clustering/all_genomes

cd $alldir

rm ./idx/*
wget -N ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xz -C ./idx/ -f taxdump.tar.gz

sqsub -q serial -r 1d -o $alldir/idx/output.updateidx ~/repo/updateidx.pl

wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.faa.tar.gz
tar -xz -C ./faa/ --strip-components 2 -f all.faa.tar.gz

wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz
tar -xz -C ./fna/ --strip-components 2 -f all.fna.tar.gz

wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.ffn.tar.gz
tar -xz -C ./ffn/ --strip-components 2 -f all.ffn.tar.gz

wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.ptt.tar.gz
tar -xz -C ./ptt/ --strip-components 2 -f all.ptt.tar.gz

wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gbk.tar.gz
tar -xz -C ./gbk/ --strip-components 2 -f all.gbk.tar.gz


cd $alldir
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/summary.txt
todate=`date +"%Y-%m-%d"`
cp summary.txt NCBI-genome-summary$todate
ln -sf /work/rec3141/gene_clustering/all_genomes/NCBI-genome-summary$todate /work/rec3141/gene_clustering/all_genomes/NCBI-genome-summary.newest

wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lproks_0.txt
wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lproks_1.txt
wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lproks_2.txt
wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lenvs.txt
wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/leuks.txt

cd $alldir/faa/
for file in *.faa
 do 
   formatdb -i $file
 done;

sqsub -q serial -r 1d -o $alldir/faa/output.accession_list ~/repo/accession_list.pl
