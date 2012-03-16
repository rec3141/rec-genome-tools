#!/bin/bash
#give it the taxon id to cluster and the taxon id cluster list
#it will run the clustering for different cutoffs and lengths
#e.g. ./cluster.sh AR_110001 taxalist-Archaea-with-plasmids
if [ ! -d "$1" ]; then mkdir $1;fi;

rm newcluster.sh
 for cutoff in -3 -5 -10 -20;
   do
   for length in .6 .7 .8;
   do
    echo "sqsub -r 10080m --mpp=7g -o ./$1/output.cluster"."$1".1e"$cutoff$length -q serial /home/rec3141/repo/cluster.pl $2 $1 1e$cutoff $length long on >>cluster_output.$1 2>&1" >> newcluster.sh
  done;  
 done;

chmod +x newcluster.sh
