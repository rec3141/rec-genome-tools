#!/bin/bash
# gets (~single-copy) genes in requested clusters and aligns them with MUSCLE
# then concatenates them and makes 2 phylogenetic trees:
# 1) with PhyML
# 2) with GBLOCKS and PhyML
#
# run in cluster directory (e.g. SRB/SF_110000/)
#
# given:
# 1) full path to list of single copy clusters (1\n2\n7\n12\n...)
# 2) full path to taxalist
# 3) full path to cluster_list
# 4) taxa id

#unset PERL5OPT

export PATH=/opt/sharcnet/jdk/current/bin:$PATH
clusters=$1
taxalist=$2
clusterlist=$3
taxid=$4

homedir=`pwd`
treedir=`pwd`/trees

if [ ! -d $treedir ]
 then mkdir $treedir
fi

if [ -f $treedir/cattree.run ]; then
	rm $treedir/cattree.run;
fi;


AlignMUSCLE() {	
	 for cluster in `cat $clusters`;
	  do
		if [ ! -f $treedir/$taxid.$cluster.aa.fasta.* ] 
			then 
			 if [ -f $treedir/$taxid.$cluster.run ] 
			  then rm $treedir/$taxid.$cluster.run
			 fi
			echo "cd $homedir; ~/repo/cluster_info.pl fastaaa 0 $clusterlist $cluster $taxalist $taxid > $treedir/$taxid.$cluster.aa" > $treedir/$taxid.$cluster.run
			echo "cd $treedir; cat $taxid.$cluster.aa | ~/program/muscle/muscle -diags -out $taxid.$cluster.aa.fasta" >> $treedir/$taxid.$cluster.run
			chmod +x $treedir/$taxid.$cluster.run
			echo "sqsub -q serial -r 2d -o $treedir/$taxid.$cluster.align.output $treedir/$taxid.$cluster.run" >> $treedir/cattree.run
		fi;
	  done;
}

AlignCLUSTALO() {
	for cluster in `cat $clusters`;
	  do
		if [ ! -f $treedir/$taxid.$cluster.clo.fasta ] 
			then 
			 if [ -f $treedir/$taxid.$cluster.run ] 
			  then rm $treedir/$taxid.$cluster.run
			 fi
			cp $treedir/$taxid.$cluster.aa.fasta $treedir/$taxid.$cluster.aa.1st
			echo "cd $treedir;  ~/bin/clustalo -v --iter=5 --outfmt=fasta --in $taxid.$cluster.aa.1st --out $taxid.$cluster.clo.fasta" > $treedir/$taxid.$cluster.run
			chmod +x $treedir/$taxid.$cluster.run
			echo "sqsub -q serial -r 2d -o $treedir/$taxid.$cluster.align.output $treedir/$taxid.$cluster.run" >> $treedir/cattree.run
		fi;
	  done;
}

#now we have alignments for each single-copy-cluster
#next we need to concatenate the aligments and convert to phylip format

FastaConcat() {
#	echo "sqsub -q serial -o $treedir/output.$taxid.concat.average -r 1d ~/repo/fastaconcat.pl $taxalist $taxid $treedir average" >> $treedir/cattree.run
#	echo "sqsub -q serial -o $treedir/output.$taxid.concat.lowest -r 1d ~/repo/fastaconcat.pl $taxalist $taxid $treedir lowest" >> $treedir/cattree.run
#	echo "sqsub -q serial -o $treedir/output.$taxid.concat.longest -r 1h ~/repo/fastaconcat.pl $taxalist $taxid $treedir longest" >> $treedir/cattree.run
	echo "~/repo/fastaconcat.pl $taxalist $taxid $treedir longest" >> $treedir/cattree.run
#	echo "sqsub -q serial -o $treedir/output.$taxid.concat.lengthdiff -r 1h ~/repo/fastaconcat.pl $taxalist $taxid $treedir lengthdiff" >> $treedir/cattree.run
}

ReadSeq() {
	ln -s $taxid.concat.average.faa $taxid.concat.faa
	cat $treedir/$taxid.concat.faa | java -jar ~/program/readseq/readseq.jar -a -f12 -p -o $treedir/$taxid.concat.phylip
}

Gblocks() {
	#next make a gblocks-filtered version
	cut -f1 -d' ' $treedir/$taxid.concat.faa > $treedir/$taxid.concat.cut.faa
	~/program/Gblocks/Gblocks $treedir/$taxid.concat.cut.faa -b5=h -e=.gb
	mv $treedir/$taxid.concat.cut.faa.gb $treedir/$taxid.gb.concat.faa 
	cat $treedir/$taxid.gb.concat.faa | java -jar ~/program/readseq/readseq.jar -a -f12 -p -o $treedir/$taxid.gb.concat.phylip
}


SerialTrees() {
	ln -sf $taxid.average.raxml.q.txt raxml.q.txt
	
	if [ ! -f $treedir/$taxid.raxml.s.concat.output]; then
# 	Multi gene alignment splitting (-f s) not implemented for the MPI-Version
		echo "cd $treedir; ~/bin/raxmlHPC-SSE3 -f s -p 31415 -b 51413 -N 3 -n $taxid.s -s $treedir/$taxid.concat.phylip -m PROTGAMMAWAG -q raxml.q.txt" > $treedir/$taxid.s.raxml
		chmod +x $treedir/$taxid.s.raxml
		echo "sqsub -q serial -r 1d --mpp=3800M -o $treedir/$taxid.raxml.s.concat.output $treedir/$taxid.s.raxml" >> $treedir/cattree.run
	fi;
	
	if [ ! -f $treedir/$taxid.phyml.output ]; then
		echo "cd $treedir; ~/bin/phyml -i $taxid.concat.phylip -d aa" > $treedir/$taxid.phyml
		chmod +x $treedir/$taxid.phyml
		echo "sqsub -q serial -r 3d --mpp=3800M -o $treedir/$taxid.phyml.output $treedir/$taxid.phyml" >> $treedir/cattree.run
	fi;

	if [ ! -f $treedir/$taxid.gb.phyml.output ]; then
		echo "cd $treedir; ~/bin/phyml -i $taxid.gb.concat.phylip -d aa --no_memory_check" > $treedir/$taxid.gb.phyml
		chmod +x $treedir/$taxid.gb.phyml
		echo "sqsub -q serial -r 3d --mpp=3800M -o $treedir/$taxid.gb.phyml.output $treedir/$taxid.gb.phyml" >> $treedir/cattree.run
	fi;
}

#then make a tree from the concatenated alignment using phyml
MPITrees() {
	ln -sf $taxid.average.raxml.q.txt raxml.q.txt

	if [ ! -f $treedir/$taxid.raxml.d.concat.output ]; then
		echo "cd $treedir; ~/bin/raxmlHPC-HYBRID-SSE3 -f d -p 31415 -b 51413 -N 100 -n $taxid.d -s $treedir/$taxid.concat.phylip -T 8 -m PROTGAMMAWAG " > $treedir/$taxid.d.raxml
		chmod +x $treedir/$taxid.d.raxml
		echo "sqsub -fmpi -r 1d -N 2 -n 16 --mpp=3800M -o $treedir/$taxid.raxml.d.concat.output $treedir/$taxid.d.raxml" >> $treedir/cattree.run
	fi;
	
	if [ ! -f $treedir/$taxid.phyml.boot.output ]; then
		echo "cd $treedir; ~/bin/phyml-mpi -i $taxid.concat.phylip -d aa -n 100 --no_memory_check" > $treedir/$taxid.boot.phyml
		chmod +x $treedir/$taxid.boot.phyml
		echo "sqsub -fmpi -r 3d -N 2 -n 16 --mpp=3800M -o $treedir/$taxid.phyml.boot.output $treedir/$taxid.boot.phyml" >> $treedir/cattree.run
	fi;
	
	if [ ! -f $treedir/$taxid.raxml.d.gb.concat.output ]; then
		echo "cd $treedir; ~/bin/raxmlHPC-HYBRID-SSE3 -f d -p 31415 -b 51413 -N 100 -n $taxid.d.gb -s $treedir/$taxid.gb.concat.phylip -T 8 -m PROTGAMMAWAG " > $treedir/$taxid.d.gb.raxml
		chmod +x $treedir/$taxid.d.gb.raxml
		echo "sqsub -fmpi -r 1d -N 2 -n 16 --mpp=3800M -o $treedir/$taxid.raxml.d.gb.concat.output $treedir/$taxid.d.gb.raxml" >> $treedir/cattree.run
	fi;
}


#AlignMUSCLE;
#AlignCLUSTALO;
#FastaConcat; 
#ReadSeq;
Gblocks;
#SerialTrees;
#MPITrees;
