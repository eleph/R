# 
#http://www.ipsj.or.jp/10jigyo/taikai/73kai/73program/data/pdf/5J-1.html
# title: A Proposal for Outlier Detection in High Dimensional Space
# Author: Zhana
# Email: znabao@ruri.waseda.jp

# pre define
m <- 100
n <- 1000
d <- 20
c <- ceiling(n/d)
s <- array(0, c((2*m),c))
t <- array(0,c((2*m),n))
dat <- read.table("ArtiData1001000.dat", sep=" ")      #data[n,m]
d <- array(0,c(m))
SI <- array(0, c(n))


#step1 calculate section number in dimension
  imax = sapply(dat, max)
  imin = sapply(dat, min)
  len = (imax - imin) * 1.001       # enlarge the scope avoid margin problem
  imax = imax + len * 0.0005
  imin = imin - len * 0.0005
  unt = len /c
  for(i in 1:m)
  {
    dmin = imin[i]
    dunt = unt[i]
    k = ceiling((dat[,i]- dmin)/dunt)
    for(j in 1:n){
            t[i,j] <- k[j]
            s[i,k[j]] <- s[i,k[j]]+1
                  } 
	ci = length(s[i,s[i,]!= 0])
	d[i] <- n/ci
  }
  
  
# Setp2 Calculate section density when points dimension change i->i+1

for (i in 1:(m-1))
for(j in 1:c)
   if(s[i,j])s[i+m,j] = sum(t[i+1,t[i,]==j])/s[i,j]
                       
for(j in 1:c)
   if(s[m,j]) s[(m+m),j] = sum(t[1,t[m,]==j])/s[m,j]
                        


 # Setp3 Grid Transform Self-information Statistic
 for(i in 1:m)
 for(j in 1:n)
 {
  val = t[i,j]
  
  if(s[i,val]<3) t[i+m,j]= 0
  else {
      
	  if(i<m){k=i+1}
	  else {k=1}
	  
	  tmp <- t[k,t[i,]== val]
	  tmax = max(tmp)
      tmin = min(tmp)
	  
	  po<- s[i+m,t[i,j]]
	  
	  if(tmax==po)
	  {
		  Pcent=0
	  }
	  else
	  {
		  Pcent <- round((po-tmin)/(tmax-po))
	  }
	if(Pcent==1)
	  {
		  if(t[i+1,j]<po)
		  {
			  tmp2 <- tmp[tmp[]<po]
			  po1<- sum(tmp2)/length(tmp2)
			  avg <- mean(abs(tmp2[]-po1))
			  pval <-abs(t[k,j]-po1)
			  if(pval!=0&& avg!=0)
			  {
				  if((1.5*avg)<pval){
					  t[i+m,j] = log2(avg/pval)
				  }
				  else
				  {
					  t[i+m,j]=0
				  }
			  }
			  else
			  {
				  t[i+m,j]=0
			  }
		  }
		else
		{
				  tmp2 <- tmp[tmp[]>po]
				  po1<- sum(tmp2)/length(tmp2)
				  avg <- mean(abs(tmp2[]-po1))
				  pval <-abs(t[k,j]-po1)
				  if(pval!=0&& avg!=0)
				  {
					  if((1.5*avg)<pval){
						  t[i+m,j] = log2(avg/pval)
					  }
					  else
					  {
						  t[i+m,j]=0
					  }
				  }
				  else
				  {
					  t[i+m,j]=0
				  }
			  }
	  }
				  
	else
	{
	  avg <- mean(abs(tmp[]-po))
	  pval <-abs(t[k,j]-po)
	 
	 
	  if(pval!=0&& avg!=0)
	  {
		  if((1.5*avg)<pval){
		  t[i+m,j] = log2(avg/pval)
		  }
		  else
		  {
			t[i+m,j]=0
		  }
	  }
	  else
	  {
		 t[i+m,j]=0
	  }
	}
	  
   }
  } 
  for(i in 1:m)
  {
    s[i,(1.5*s[i,])>d[i]] =d[i]
    for(j in 1:n) 
      {
		  
		  t[i,j] = log2(s[i,t[i,j]]/d[i])
		  
	  }
  }
  
  
  
  sapply(1:n, function(j)SI[j] <<- sum(t[,j ]))
  
  shreshod <- -60
  outlier <- SI[SI<(shreshod)]
  position <- which(SI<(shreshod))
  output <- cbind(outlier,position)

  x <- c(1:n)
  plot(x, SI)

  