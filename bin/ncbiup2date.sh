#!/bin/bash

#run in SCREEN not on SQSUB because compute nodes can't access internet

UpdateTaxonomy() {
echo "Updating Taxonomy..."
wget -N ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
rm ./idx/*
tar -xz -C ./idx/ -f taxdump.tar.gz
echo "  done Updating Taxonomy."

echo "Submitting BioPerl taxonomy update to cluster..."
sqsub -q serial -r 1d -o $alldir/idx/output.updateidx ~/repo/updateidx.pl
echo "  submitted."
}

DownloadComplete() {
echo "Downloading Bacterial genome information..."
wget -N http://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/summary.txt
todate=`date +"%Y-%m-%d"`
cp summary.txt NCBI-genome-summary$todate
ln -sf /work/rec3141/gene_clustering/all_genomes/NCBI-genome-summary$todate /work/rec3141/gene_clustering/all_genomes/NCBI-genome-summary.newest

#these files are deprecated
#wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lproks_0.txt
#wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lproks_1.txt
#wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lproks_2.txt
#wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/lenvs.txt
#wget -N ftp://ftp.ncbi.nih.gov/genomes/genomeprj/leuks.txt

echo "Downloading microbial genome 'faa' files"
wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.faa.tar.gz
tar -xz -C ./faa/ --strip-components 2 -f all.faa.tar.gz

echo "Downloading microbial genome 'fna' files"
wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz
tar -xz -C ./fna/ --strip-components 2 -f all.fna.tar.gz

echo "Downloading microbial genome 'ffn' files"
wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.ffn.tar.gz
tar -xz -C ./ffn/ --strip-components 2 -f all.ffn.tar.gz

echo "Downloading microbial genome 'ptt' files"
wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.ptt.tar.gz
tar -xz -C ./ptt/ --strip-components 2 -f all.ptt.tar.gz

#wget -N ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gbk.tar.gz
#tar -xz -C ./gbk/ --strip-components 2 -f all.gbk.tar.gz
}

DownloadDrafts() {
cd $alldir
echo "Downloading draft microbial genome files"
lftp -e 'open ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/; 
  nlist */*contig.faa*.tgz > contig.faa.draftlist; 
  nlist */*scaffold.faa*.tgz > scaffold.faa.draftlist; 
  nlist */*contig.fna*.tgz > contig.fna.draftlist; 
  nlist */*scaffold.fna*.tgz > scaffold.fna.draftlist; 
  nlist */*contig.ffn*.tgz > contig.ffn.draftlist; 
  nlist */*scaffold.ffn*.tgz > scaffold.ffn.draftlist; 
  nlist */*contig.ptt*.tgz > contig.ptt.draftlist; 
  nlist */*scaffold.ptt*.tgz > scaffold.ptt.draftlist; 
bye'

rm $alldir/contig/contig.run
rm $alldir/scaffold/scaffold.run

for tgzlist in *.draftlist
  do 
   cd $alldir
   echo "Processing $tgzlist"

   ftype=`echo $tgzlist | cut -f1 -d'.'`
   ext=`echo $tgzlist | cut -f2 -d'.'`
   sed -e 's|.*|ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/&|' $tgzlist > $tgzlist.wget

   cd $alldir/$ftype/$ext
   rm -f $alldir/$ftype/$ext/*.$ext
   wget -N -r -l inf -np -nd -i $alldir/$tgzlist.wget

   for tgz in `cat $alldir/$tgzlist`
    do 
	 filename=`echo $tgz | cut -f2 -d'/'`
	 refseq=`echo $filename | cut -f1 -d'.'`
	 echo "tar xOzf $alldir/$ftype/$ext/$filename >> $alldir/$ftype/$ext/$refseq.$ext" >> $alldir/$ftype/$ftype.run
    done

  done

   
  
}

MakeBLASTdb() {
	sqsub -q serial -o output.makeblastdb -r 7d find $alldir -name "*.faa" -exec ~/program/bin/makeblastdb -dbtype prot -in {} \;
}


UpdateAccessionList() {
	sqsub -q serial -r 1d -o $alldir/faa/output.accession_index ~/repo/accession_index.pl
	sqsub -q serial -r 1d -o $alldir/faa/output.accession_list ~/repo/accession_list.pl
}

#############
alldir=/work/rec3141/gene_clustering/all_genomes
cd $alldir

#UpdateTaxonomy;
#DownloadComplete;
#DownloadDrafts;
#MakeBLASTdb;
#UpdateAccessionList;

exit 0;
