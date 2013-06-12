#!/bin/bash
#give it the taxon id to cluster and the taxon id cluster list
#it will run the clustering for different cutoffs and lengths
#e.g. ./cluster.sh AR_110001 taxalist-Archaea-with-plasmids
if [ ! -d "$1" ]; then mkdir $1;fi;

if [ -f "newcluster.sh" ]; then rm newcluster.sh; fi;

 for cutoff in -3 -5 -10 -20 -50;
   do
   for length in .6 .7 .8;
   do
    echo "sqsub -r 10080m --mpp=15g -o ./$1/output.cluster"."$1".1e"$cutoff$length.link -q serial /home/rec3141/repo/cluster.pl $2 $1 1e$cutoff $length seed on link >>cluster_output.$1 2>&1" >> newcluster.sh
  done;  
 done;

for inflation in 1.4 2 4 6;
  do
  for cutoff in -3 -10 -20;
  do
  for length in .6 .8;
  do
    echo "sqsub -r 10080m --mpp=31g -o ./$1/output.cluster"."$1".1e"$cutoff$length.mcl_$inflation -q serial /home/rec3141/repo/cluster.pl $2 $1 1e$cutoff $length none on mcl $inflation >>cluster_output.$1 2>&1" >> newcluster.sh
  done;
  done;
done;

chmod +x newcluster.sh


# [rec3141@bul131 SF_110000]$ mcl --version
# mcl 12-068

#[rec3141@bul131 SF_110000]$ mcl -z
# [mcl] cell size: 8
# [mcl] cell contents: int and float
# [mcl] largest index allowed: 2147483647
# [mcl] smallest index allowed: 0
# Prune number                               10000        [-P n]
# Selection number                            1100        [-S n]
# Recovery number                             1400        [-R n]
# Recovery percentage                           90        [-pct n]
# warn-pct                                      10        [-warn-pct k]
# warn-factor                                 1000        [-warn-factor k]
# dumpstem                                                [-dump-stem str]
# Initial loop length                            0        [-l n]
# Main loop length                           10000        [-L n]
# Initial inflation                              2.0      [-i f]
# Main inflation                                 2.0      [-I f]

# [rec3141@bul131 SF_110000]$ mcl --show-schemes
#              Pruning      Selection       Recovery  Recover percentage
# Scheme 1        3000            400            500             90
# Scheme 2        4000            500            600             90
# Scheme 3        5000            600            700             90
# Scheme 4        6000            700            800             90
# Scheme 5        7000            800            900             90
# Scheme 6       10000           1100           1400             90  <--- default
# Scheme 7       10000           1200           1600             90
