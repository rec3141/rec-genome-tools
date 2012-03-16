#!/bin/bash
prokdir="/work/rec3141/gene_clustering/prok"
refdir="/work/rec3141/gene_clustering/all_genomes/faa"
scafdir="/work/rec3141/gene_clustering/all_genomes/scaffold/faa"
contdir="/work/rec3141/gene_clustering/all_genomes/contig/faa"
metadir="/work/rec3141/gene_clustering/all_genomes/metagenomes/faa"

if [ -f $refdir/$1.faa ];
then dir1=$refdir
elif [ -f $scafdir/$1.faa ];
then dir1=$scafdir
elif [ -f $contdir/$1.faa ];
then dir1=$contdir
elif [ -f $metadir/$1.faa ];
then dir1=$metadir
else 
echo "$1.faa not found"; exit;
fi

if [ -f $refdir/$2.faa ];
then dir2=$refdir
elif [ -f $scafdir/$2.faa ]; 
then dir2=$scafdir
elif [ -f $contdir/$2.faa ];
then dir2=$contdir
elif [ -f $metadir/$2.faa ];
then dir2=$metadir
else 
echo "$2.faa not found"; exit;
fi

if [ ! -d $prokdir/$1 ]; then mkdir $prokdir/$1; fi

if [ ! -f $prokdir/$1/$1.$2.hits ];
 then 
  echo "blasting $1 $2"
  ~/bin/blast-2.2.21/bin/blastall -p blastp -m7 -d $dir1/$1.faa -i $dir2/$2.faa -o $prokdir/$1/$1.$2.xml.blastout >> blast_output 2>&1;
  echo "extracting $1 $2"
  /home/rec3141/repo/hit_extractor.pl $1 $2 0.001 >> extractor_output 2>&1 ;
  rm $prokdir/$1/$1.$2.xml.blastout
fi
