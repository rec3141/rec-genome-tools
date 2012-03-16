#!/bin/bash
#include taxalist as only parameter
if [ ! "$#" -eq 1 ]
 then echo "provide taxalist"
 exit 1
fi

if [ ! -f $1 ]
 then echo "could not find taxalist $1"
 exit 1
fi


for file in `find ./ -name "*.cluster_list" | grep -v "gm"`;
	do
	fbase=`basename $file`;
	fdir=`dirname $file | xargs basename`;
     if [ ! -f output.gm.$fbase ]
    then echo sqsub -q serial -r 10080m --mpp=14g -o $fdir/output.gm.$fbase /home/rec3141/repo/genematrix.pl $1 $file $fdir;
  fi
done
