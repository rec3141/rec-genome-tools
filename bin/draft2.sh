alldir=/work/rec3141/gene_clustering/all_genomes/
cd $alldir

lftp -e 'open ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/; 
  nlist */*contig.faa*.tgz > contig.faa.list; 
  nlist */*scaffold.faa*.tgz > scaffold.faa.list; 
  nlist */*contig.fna*.tgz > contig.fna.list; 
  nlist */*scaffold.fna*.tgz > scaffold.fna.list; 
  nlist */*contig.ffn*.tgz > contig.ffn.list; 
  nlist */*scaffold.ffn*.tgz > scaffold.ffn.list; 
  nlist */*contig.ptt*.tgz > contig.ptt.list; 
  nlist */*scaffold.ptt*.tgz > scaffold.ptt.list; 
#  nlist */*contig.gbk*.tgz > contig.gbk.list; 
#  nlist */*scaffold.gbk*.tgz > scaffold.gbk.list; 
bye'

for file in *.list
  do 
   ftype=`cut -f1 -d'.' $file`
   ext=`cut -f2 -d'.' $file`

  for tgz in `echo $file`
   do 
     wget -N -r -l inf -np -P $ftype/$ext/ -nd ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/$tgz
   done
  done
