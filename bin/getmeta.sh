#/bin/bash
# input 1: blastlist of project accessions; e.g. 4567899.3 (mg-rast), 2001200000 (jgi-m), 637000001 (jgi-w)
# input 2: database ('mg-rast', 'img-m', 'img-w')
# input 3: required for 'mg-rast': metadata mapping file; metadata.ids.csv

metafile=$1
database=$2
metadata=$3

for id in `cat $metafile | cut -f1`;
 do
  if [ "$database" == "mg-rast" ]; then
    idd=`echo $id | tr '.' '_'`
#    prj=`grep -m1 -P "^\Q$id\E" $metadata | cut -f3`
    prj=`grep -m1 -P "^\Q$id\E" $metadata | cut -f2 -d' '`
    echo "$database: $id $prj"
    if [ ! -n "$prj" ]; then continue; fi

    wget -N -P processed/ ftp://ftp.metagenomics.anl.gov/metagenomes/$prj/$id.processed.tar.gz

    if [ ! -n "processed/$id.processed.tar.gz" ]; then continue; fi;

    filedate=$(stat --format '%Z' processed/$id.processed.tar.gz)
    nowdate=$(date +%s)
    if (( $nowdate - $filedate <6000 )); then
      tar xvzf processed/$id.processed.tar.gz $id/processed/299.screen.passed.fna.gz -O | gunzip > fna/MG_"$idd".fna
      tar xvzf processed/$id.processed.tar.gz $id/processed/350.genecalling.coding.fna.gz -O | gunzip > fna/MG_"$idd".genes.fna
      tar xvzf processed/$id.processed.tar.gz $id/processed/350.genecalling.coding.faa.gz -O | gunzip > faa/MG_"$idd".genes.faa
    fi

  elif [ "$database" == "img-m" ];
   then
    echo "$database: $id"
    idd=`echo $id | tr '.' '_'`
    wget -N -P processed/ ftp://ftp.jgi-psf.org/pub/IMG/img_mer_v340/$id.tar.gz

    if [ ! -n "processed/$id.tar.gz" ]; then continue; fi;
    filedate=$(stat --format '%Z' processed/$id.tar.gz)
    nowdate=$(date +%s)
    if (( $nowdate - $filedate <6000 )); then
      tar xvzf processed/$id.tar.gz $id.fna -O > fna/IMGM_"$idd".fna
      tar xvzf processed/$id.tar.gz $id.genes.fna -O > fna/IMGM_"$idd".genes.fna
      tar xvzf processed/$id.tar.gz $id.genes.faa -O > faa/IMGM_"$idd".genes.faa
    fi

  elif [ "$database" == "img-w" ];
   then
    echo "$database: $id"
    idd=`echo $id | tr '.' '_'`

    wget -N -P processed/ ftp://ftp.jgi-psf.org/pub/IMG/img_w_v350/$id.tar.gz

    if [ ! -n "processed/$id.tar.gz" ]; then continue; fi;
    filedate=$(stat --format '%Z' processed/$id.tar.gz)
    nowdate=$(date +%s)
    if (( $nowdate - $filedate <6000 )); then
      tar xvzf processed/$id.tar.gz $id.fna -O > fna/IMGW_"$idd".fna
      tar xvzf processed/$id.tar.gz $id.genes.fna -O > fna/IMGM_"$idd".genes.fna
      tar xvzf processed/$id.tar.gz $id.genes.faa -O > faa/IMGW_"$idd".genes.faa
    fi;
  fi;
done;
