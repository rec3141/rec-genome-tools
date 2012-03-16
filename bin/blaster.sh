#!/bin/bash
for file1 in *.faa;
 do f1s=`echo $file1 | cut -f1 -d\.`;
 for file2 in *.faa; 
  do f2s=`echo $file2 | cut -f1 -d\.`; 
   if [ ! -f $f1s"."$f2s.hits ];
   then sqsub -q serial -r 1h -o $1"."$2.jobout ./signifier.sh $f1s $f2s $file1 $file2 >> submission_record 2>&1;
   fi; 
 done; 
done