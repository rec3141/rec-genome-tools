#/bin/bash
#give 2 inputs: 1) blastlist of whole genomes, 2) blastlist of metagenomes
rm base.run
wgsdir='/work/rec3141/gene_clustering/all_genomes/fna'
metadir='/work/rec3141/gene_clustering/all_genomes/metagenomes/fna'

for f1 in `cat $1`
 do file1="$f1.fna"
  for f2 in `cat $2`
   do file2="$f2.fna"
   echo "~/program/bwa-0.6.1/bwa bwasw $wgsdir/$file1 $metadir/$file2 > $f1.$f2.sam" >> base.run
   done;
 done;

rm dobwa.run
split -l 50 -a 5 -d base.run bwa
chmod +x bwa*
for bwa in bwa*
 do echo "sqsub -q serial -o output.$bwa -r 7d ./$bwa" >> dobwa.run
done;
chmod +x dobwa.run

