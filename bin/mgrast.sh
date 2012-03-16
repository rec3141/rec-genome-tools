#/bin/bash
# input 1: blastlist of mg-rast projects
# input 2: metadata csv file
for id in `cat $1`
 do
  $id2=`echo $id | sed 's/\./\_/'
  prj=`grep -m1 -P "^\Q$id2\E" $2 | cut -f2`
  wget -N -P raw/ ftp://ftp.metagenomics.anl.gov/metagenomes/$prj/$id.raw.tar.gz
  wget -N -P processed/ ftp://ftp.metagenomics.anl.gov/metagenomes/$prj/$id.processed.tar.gz
  tar xvzf processed/$id.processed.tar.gz $id/processed/*genecalling.coding.fna.gz -O | gunzip > fna/MG_$id2.fna
  tar xvzf processed/$project.processed.tar.gz $project/processed/*genecalling.coding.faa.gz -O | gunzip > faa/MG_$id2.faa
 done;