#for file in `find $1 -name "*.hits"`; 
for file in `cat $1`;
 do 
  basefile=`basename $file`;
  prokdir=`echo $basefile | cut -f1 -d'.'`;
#  ls -l $file
#  ls -l  /work/rec3141/gene_clustering/prok/$prokdir/$basefile
  if [[ $file -nt /work/rec3141/gene_clustering/prok/$prokdir/$basefile ]]; then
   mv -v $file /work/rec3141/gene_clustering/prok/$prokdir/$basefile;
  else
   rm -v $file;
 fi;
done;

