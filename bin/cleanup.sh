#!/bin/bash
list=`find -L ./$1 -name "*xml*"`;
for f1 in $list;
 do 
 mydir=`dirname $f1`
 myfile=`basename $f1`
 
 f1s=`echo $myfile | cut -f1 -d\.`;
 f2s=`echo $myfile | cut -f2 -d\.`;

 if [ -f $mydir/$f1s.$f2s.hits ]
   then echo rm $mydir/$f1s.$f2s.hits; 
 fi
 if [ -f $mydir/$myfile ]
   then echo rm $mydir/$myfile
 fi

done
