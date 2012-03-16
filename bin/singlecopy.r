fl <- list.files(path=".",pattern="*cluster_list$")

for (myclust in fl) {
  myfile <- paste(myclust,".r.csv",sep="")
  print(myfile)
  a<-read.csv(myfile,sep="\t",head=T,row.names=1)
  b<-a[,3:(ncol(a)-4)]
  c<-which(rowSums(b>0)==rowSums(b) & rowSums(b)==ncol(b))
  print(length(c))
  write.table(c,file=paste(myclust,".single-copy.csv",sep=""),row.names=F,col.names=F)
}