#!/bin/bash
# given:
# 1) full path to list of single copy clusters (1\n2\n7\n12\n...)
# 2) full path to taxalist
# 3) full path to cluster_list
# 4) taxa id
export PATH=/opt/sharcnet/jdk/current/bin:$PATH
clusters=$1
taxalist=$2
clusterlist=$3
taxid=$4
treedir=`pwd`/trees

if [ ! -d $treedir ]
 then mkdir $treedir
fi

#myfile=`basename $clusterlist`
myfile=$taxid

if [ -f $treedir/$myfile.phyml ]
 then
    rm $treedir/$myfile.phyml
    touch $treedir/$myfile.phyml
fi

 for cluster in `cat $clusters`;
  do
	if [ ! -f $treedir/$myfile.$cluster.aa.fasta ] 
	then 
	 if [ -f $treedir/$myfile.$cluster.run ] 
	  then rm $treedir/$myfile.$cluster.run
         fi
	echo "~/repo/retrieve_sequences.pl fastaaa 0 $cluster $clusterlist $taxalist $taxid > $treedir/$myfile.$cluster.aa" > $treedir/$myfile.$cluster.run
#	echo "~/repo/retrieve_sequences.pl $clusterlist $cluster fastant $taxalist $taxid > $treedir/$myfile.$cluster.nt" >> $treedir/$myfile.$cluster.run
	echo "cd $treedir; cat $myfile.$cluster.aa | ~/program/muscle/muscle -diags -stable -out $myfile.$cluster.aa.fasta" >> $treedir/$myfile.$cluster.run
	chmod +x $treedir/$myfile.$cluster.run
	sqsub -q serial -r 7d --mpp=7g -o $treedir/$myfile.$cluster.concat.output $treedir/$myfile.$cluster.run
	fi;
    done;

#wait for them all to finish
while [ `ls $treedir | grep -c run` -ne `ls $treedir | grep -c fasta` ]
do sleep 60
done

#now we have alignments for each single-copy-cluster
#next we need to concatenate the aligments and convert to phylip format
~/repo/fastaconcat.pl $treedir $taxid $taxalist

cat $treedir/$myfile.concat.faa | java -jar ~/program/readseq/readseq.jar -a -f12 -p -o $treedir/$myfile.concat.phylip
#then make a tree from the concatenated alignment using phyml
echo "cd $treedir; ~/program/phyml_3.0.1_linux64 -i $myfile.concat.phylip -d aa" > $treedir/$myfile.phyml
chmod +x $treedir/$myfile.phyml
sqsub -q serial -r 7d -o $treedir/$myfile.concat.output --mpp=31g $treedir/$myfile.phyml

#next make a gblocks-filtered version
~/program/Gblocks/Gblocks $treedir/$myfile.concat.faa -b5=h -e=.gb
mv $treedir/$myfile.concat.faa.gb $treedir/$myfile.gb.concat.faa 
cat $treedir/$myfile.gb.concat.faa | java -jar ~/program/readseq/readseq.jar -a -f12 -p -o $treedir/$myfile.gb.concat.phylip
#then make a tree from the concatenated alignment using phyml
echo "cd $treedir; ~/program/phyml_3.0.1_linux64 -i $myfile.gb.concat.phylip -d aa" > $treedir/$myfile.gb.phyml
chmod +x $treedir/$myfile.gb.phyml
sqsub -q serial -r 7d -o $treedir/$myfile.gb.concat.output --mpp=31g $treedir/$myfile.gb.phyml


