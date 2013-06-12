for file in *.sam
 do 
  f1=`basename $file .sam`
  samtools view -b -F 0x04 -S $file > $f1.bam
 done
