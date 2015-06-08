#Title: Two Phases Outlier Detection in Different Subspaces
#PhD Workshop PIKM at CIKM 2014
#http://dl.acm.org/citation.cfm?id=2668046#references
# Author: Zhana
# Email: znabao@ruri.waseda.jp

library(ROCR)
dm = 100
pn = 500
K1=10
K2 = 4
#K2 < K1 * 2/3
lp=5  #映射维数
#k = 5+ceiling(runif(1)*5)
nnd <-array(0, c((K1+1),pn))   # 邻居集
knn <-array(0, c((K1+1),pn))
denk <-array(0, c((K1+1)))

den <-array(0, c(pn,dm))
dev <-array(0, c(pn,dm))

rdm <- c(1:dm)

dat<- read.table("ArtiData100500.dat", sep=" ")
dat<- as.matrix(dat)


tim1<- proc.time()

for(j in 1:dm)
{
 
for(i in 1:pn)
{
  #diff <- abs(scale(dat[,j],dat[i,j],FALSE))
  dif <- abs(dat[,j]-dat[i,j])    
  order.dist <- order(dif) #dtmp
  orddist <- order.dist[1:(K1+1)]
  nnd[,i] <- orddist
  orddist <- orddist[-1]
  tmp <- dif[orddist]
  
  den[i,j] <- (mean(tmp)+ tmp[K1])/2
   
}
# 映射 密度
 den[,j] <-mean(den[,j])/den[,j] 
 #den[den[,j]>1,j]=1
 
 
 rdm <- c(1:dm)
 rdm <- rdm[-j]
 rd <- sample(rdm,size=lp)     # 维数大于5时
 #id <-which(rdm==rd)
 #rdm <- rdm[-id]
 

 for(v in 1:pn)
	{
   if(den[v,j]<0.25)       # too sparse to be deviated
     {dev[v,j]=0
      next}
   for(u in rd)
   {   
    knn <- dat[nnd[,v],u]
    for(x in 1:(K1+1))
    {
     #diff <- abs(scale(knn[],knn[x],FALSE))
     dif <- abs(knn-knn[x])
     #order.dist <- order(diff) #dtmp
     #orddist <- order.dist[2:(K2+1)]
     #tmp <- diff[orddist]
     #denk[x] <- (mean(tmp)+ tmp[K2])/2
     tmp <- sort(dif)
     denk[x] <- (mean(tmp[2:(K2+1)])+ tmp[(K2+1)])/2
    }
    denk[1] <-mean(denk[])/denk[1]
    print(denk)
    #if(denk[1]>1)denk[1]<-1
    dev[v,j] <- dev[v,j] + log2(denk[1])
   } #end u
   dev[,j] <- dev[,j] 
  
  }  #end v 

}  # end j
 
SI <- sapply(1:pn, function(i)(-1)*(sum(log2(den[i,]))+sum(dev[i,]))/dm)



tim2<- proc.time()
tim <- tim2-tim1
print(tim)


SIod <- rev(order(SI))
outlier <- SI[SIod]
output <- cbind(outlier,SIod)
print(output[1:20,])


st<-paste("ArtiRst",dm,pn,".dat", sep = "")
z<-pn-10
rst<-output
rst[,2]=0
rst[output[,2]>z,2]=1
write.table(rst,file=st, row.names=F,col.names =F)






