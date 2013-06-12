fl <- list.files(path=".",pattern="^[A-Z]{2}_.*cluster_list.r.csv$")

for (myclust in fl) {
  myfile <- paste(myclust,".1k.csv",sep="")
  print(myfile)
  a<-read.csv(myfile,sep="\t",head=T,row.names=1)
  b<-a[,3:(ncol(a)-4)]
  ncol(b)
#  c<-which(rowSums(b>0)==rowSums(b) & rowSums(b)==ncol(b))
#  c<-which(rowSums(b>0)==rowSums(b))
#  c<-which(rowSums(b>0)==rowSums(b) & rowSums(b)>=0.95*ncol(b))
#  c<-which(rowSums(b>0)>0.9*rowSums(b) & rowSums(b)==ncol(b))
  c<-which(rowSums(b>0)>0.97*rowSums(b) & rowSums(b)>.97*ncol(b))
  d<-cbind(c,rowSums(b[c,]))
  print(dim(d))
  print(d)
  write.table(d,file=paste(myclust,".single-copy.csv",sep=""),row.names=F,col.names=F)
}
