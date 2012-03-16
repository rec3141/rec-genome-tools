keys <- read.table("key-bacillus",sep="\t",row.names=1)

fl <- list.files(path=".",pattern="^BL.+cluster_list$")
for (myclust in fl) {
  myfile <- paste(myclust,".r.csv",sep="")
  print(myfile)
  a<-read.csv(myfile,sep="\t",head=T,row.names=1)
  b<-a[,3:(ncol(a)-4)]
  mycols <- unique(unlist(lapply(keys[,2],function(x) grep(x,colnames(b)))))
  b<-b[,mycols]
  c<-which(rowSums(b>0)==rowSums(b) & rowSums(b)==ncol(b))
  write.table(c[1:100],file=paste(myclust,".single-copy.csv",sep=""),row.names=F,col.names=F)
}
