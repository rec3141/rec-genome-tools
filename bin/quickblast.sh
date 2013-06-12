#/usr/bin/bash
# file of lines to execute on sharcnet system

CheckQueue() {  

	queued=`sqjobs -q | grep -c 'quickblast'`
  	running=`sqjobs -r | grep -c 'quickblast'`

  	if [ $running -gt 0 ]; then
   		if [ $queued -lt 10 ]; then
   			sqsub -q serial -r 1d -o output.quickblast -j quickblast ~/repo/quickblast.sh $file $blastlist 0 $queue;
   		else
   			sleep 15;
   		fi;
   	else
   		sleep 15;
  	fi;
}


ExecLine() {
	workline=`perl -e '$if=shift(@ARGV);open(IF,"<$if") or die;@lines=<IF>;print shift(@lines);close(IF);open(OF,">$if"); print OF @lines; close(OF);' $file.$syst`
 	CMD=`echo $workline | tr ';' '\n'`
 	echo "begin on $syst: $CMD"
 	eval "$CMD"
 	if [ "$?" -ne "0" ]; then 
	 	echo "error on $syst: $CMD"
	 	echo $workline >> errors.$file
		#exit $?
  	else
   		echo "complete on $syst: $CMD"
  	fi;
}

GetBlastList() {
    if [ -f "$file.wait" ]; then
        sleep 10;
    else 
        touch $file.wait
	perlCMD='use File::Copy;$if=shift(@ARGV);$syst=shift(@ARGV);$queue=shift(@ARGV);open(IF,"<",$if) or die $!;@lines=<IF>;close(IF) or die $!;open(OF,">","$if.$$.tmp");print OF @lines[$queue..$#lines];close(OF) or die $!;move("$if.$$.tmp","$if");open(OF,">>","$if.$syst"); print OF @lines[0..$queue]; close(OF) or die $!;'
	perl -e "$perlCMD" $file $syst $queue
	 if [ "$?" -ne "0" ]; then 
		echo "error on $syst: perl -e \"$perlCMD\" $file $syst $queue"
	 	echo $workline >> errors.$file
	  else
   		echo "complete on $syst: perl -e \"$perlCMD\" $file $syst $queue"
  	 fi;
	rm $file.wait
    fi;
}


RunBlasts() {

	if [ $master -ge 1 ]; then
		while [ -s $file ];
	 		do
	 			CheckQueue;
 	 		done
	 	 if [ $master -gt 1 ]; then
			CleanUp;
		 fi;
  	else
		GetBlastList; #make working blastlist for each node with $queue per node
    	while [ -s $file.$syst ];
      		do
         		ExecLine;
      		done
	fi;
}

CleanUp() {
#	find -L ./ -empty -name "*.hits" -delete
#	find -L ./ -name "*.hits.tmp" -print -delete
	~/repo/blastlist_process.pl $blastlist
	cat N*.run > $file.tmp
	rm *.run
	mv $file.tmp $file
}

###############
file=$1
blastlist=$2
master=$3 #0 = spawned; 1 = spawn; 2 = do cleanup
queue=$4 #size of queue
syst=`uname -n`

#wait random amount of time between starts
number=$RANDOM 
let "number %= 12"
sleep $number

if [ $master -gt 0 ]; then
	sed -i -e '/^$/d' $file;
	sort -u $file > $file.tmp;
	mv $file.tmp $file;
fi;

if [ $master -gt 1 ]; then
	CleanUp;
fi;

while [ -s $file ];
 do
  RunBlasts; #loop de loop
 done


