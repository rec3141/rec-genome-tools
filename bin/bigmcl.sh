#!/bin/bash

function gethits () {
for file in `find ./ | grep hits`
do
cat $file | cut -f1,2,3, -d" "
done;
}

dr=`pwd | cut -f5 -d\/`

gethits | mcl - --abc --abc-neg-log -o /scratch/rec3141/bigmcl.$dr.mcl
