#!/bin/bash
for clusterfile in `ls ./AR_110001/*.cluster_list`
do
  file=`echo $clusterfile | cut -f3 -d\/`
  if [ ! -f /scratch/rec3141/$file.mcl.4.output ]
  then
  rm /scratch/rec3141/$file.mcl.4.run
  touch /scratch/rec3141/$file.mcl.4.run
  chmod +x /scratch/rec3141/$file.mcl.4.run

  numclust=`cat $clusterfile | grep -c -P "\d \D"`
  for cluster in `seq $numclust -1 1`
  do
  if [ ! -f /scratch/rec3141/$file.$cluster.mcl ]
  then echo "sed '$cluster""q;d' $clusterfile.evalues.output | tr '|' '\n' | mcl - --abc --abc-neg-log -o /scratch/rec3141/$file.$cluster.mcl" >> /scratch/rec3141/$file.mcl.4.run
#echo "~/repo/retrieve_sequences.pl $clusterfile $cluster evalue ./taxalist-Archaea-with-plasmids | mcl - --abc --abc-neg-log -o /scratch/rec3141/$file.$cluster.mcl" >> /scratch/rec3141/$file.mcl.3.run
  fi;
  done;
  echo sqsub -q serial -r 10080m -o /scratch/rec3141/$file.mcl.4.output --mpp=8g /scratch/rec3141/$file.mcl.4.run
  fi;
done;
