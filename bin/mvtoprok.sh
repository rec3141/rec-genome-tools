for file in `find $1 -name "*.hits"`; 
 do 
  basefile=`basename $file`;
  prokdir=`echo $basefile | cut -f1 -d'.'`;
  mv -v $file /work/rec3141/gene_clustering/prok/$prokdir/$basefile;
done;

