#/bin/bash
#give 2 inputs: 1) blastlist of whole genomes, 2) blastlist of metagenomes
#output to > bwa.run
dir1='/work/rec3141/gene_clustering/all_genomes/contig/fna'
dir2='/work/rec3141/gene_clustering/all_genomes/scaffold/fna'
dir3='/work/rec3141/gene_clustering/all_genomes/fna'
#prokdir='/work/rec3141/gene_clustering/prok'
prokdir='/scratch/rec3141'
metadir='/work/rec3141/gene_clustering/all_genomes/metagenomes/fna'

blastwgs=$1
blastmeta=$2

for wgsdir in $dir1 $dir2 $dir3
 do
 for f1 in `cat $blastwgs`
   do
   if [ ! -e "$wgsdir/$f1.fna" ]; then continue; fi
   for f2 in `cat $blastmeta`
     do
     if [ -e "$metadir/$f2.genes.fna" ]; then genes="genes"; file2="$f2.genes.fna"
     elif [ -e "$metadir/$f2.fna" ]; then genes="whole"; file2="$f2.fna"
     else continue;
     fi

     if [ -e "$prokdir/$f1.$f2.$genes.sam" ]; then  continue;
     else echo "~/program/bwa-0.6.1/bwa bwasw $wgsdir/$f1.fna $metadir/$file2 > $prokdir/$f1.$f2.$genes.sam"
     fi

     done;
    done;
 done;

