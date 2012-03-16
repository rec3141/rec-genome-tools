for mydir in `ls -d /work/rec3141/gene_clustering/prok/*/`
do
cd $mydir
d1=`echo $mydir | cut -f6 -d\/`
echo "tar czf /scratch/rec3141/$d1.hits.tgz *.hits" > ./targz.run
chmod +x ./targz.run
sqsub -q serial -r 1h -o %j ./targz.run
done;

