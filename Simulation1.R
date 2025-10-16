
#Code for simulation 1
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
Rcpp::sourceCpp("Standard_LC.cpp")

L <- 4

pl<-c(0.4,0.3,0.2,0.1)
c1=2
c2=3
c3=4

p1<-matrix(c(0.8,0.9,0.3,0.1,0.2,0.1,0.7,0.9),nrow=c1,ncol=L,byrow=T)
p2<-matrix(c(0.7,0.8,0.2,0.3,0.3,0.2,0.8,0.7),nrow=c1,ncol=L,byrow=T)
p3<-matrix(c(0.7,0.8,0.1,0.15,0.1,0.05,0.7,0.65,0.2,0.15,0.2,0.2),nrow=c2,ncol=L,byrow=T)
p4<-matrix(c(0.8,0.1,0.7,0.2,0.05,0.7,0.1,0.7,0.15,0.2,0.2,0.1),nrow=c2,ncol=L,byrow=T)
p5<-matrix(c(0.65,0.1,0.6,0.2,0.1,0.65,0.1,0.6,0.1,0.1,0.2,0.1,0.15,0.15,0.1,0.1),nrow=c3,ncol=L,byrow=T)
p6<-matrix(c(0.75,0.2,0.7,0.1,0.05,0.6,0.1,0.7,0.05,0.05,0.1,0.15,0.15,0.15,0.1,0.05),nrow=c3,ncol=L,byrow=T)


S <- 100
for(s in 1:S){

n=500
  
  data <- matrix(nrow=n,ncol=1) 
  data[,1] <- seq(1:n)
  
  
  data2<-NULL
  samples2 <- NULL
      samples <- sample(1:L,prob=pl,size=n,replace=TRUE) 
      for(i in (1:n)){
        data2p = cbind(sample(0:(c1-1),prob=p1[,samples[i]],size=1),sample(0:(c1-1),prob=p2[,samples[i]],size=1),sample(0:(c2-1),prob=p3[,samples[i]],size=1),sample(0:(c2-1),prob=p4[,samples[i]],size=1),sample(0:(c3-1),prob=p5[,samples[i]],size=1),sample(0:(c3-1),prob=p6[,samples[i]],size=1))
        data2=rbind(data2,data2p)
      } 
      samples2<-c(samples2,samples)  
          
    
  data<-cbind(data,data2)  
  datac<-cbind(data,samples2)

s<-data
s<-as.matrix(s)

y<-list(s[1,])
for(i in (2:dim(s)[1])){
  y[[i]]<-s[i,]
}

main2(y)


}
