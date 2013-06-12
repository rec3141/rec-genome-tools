#!/bin/bash

alldir='/work/rec3141/gene_clustering/all_genomes'
cd $alldir
for ftype in contig scaffold
  do
  if [ ! -d $alldir/$ftype ]
    then mkdir $alldir/$ftype
  fi
  cd $alldir/$ftype
  for ext in fna faa ptt ffn gbk
    do
    if [ ! -d $ext ]
      then mkdir $ext
    fi
    wget -N -r -l inf -np -P $ext/ -nd -A $ftype.$ext.tgz,$ftype.$ext.1.tgz,$ftype.$ext.2.tgz,$ftype.$ext.3.tgz ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/
    cd $alldir/$ftype/$ext
    for file in *.tgz
      do f1=`echo $file | cut -f1 -d'.'`
      if [ -f $f1.$ext ]
	then rm $f1.$ext
      fi 
      tar xOzf $file >> $f1.$ext
      if [ $ext == 'faa' ]
	then formatdb -i $ext/$f1.$ext
      fi
    done;
  cd $alldir/$ftype
  done;
cd $alldir
done;
