

#Code for simulation 2
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
Rcpp::sourceCpp("Multilevel_LC.cpp")

L <- 4

H0<-3
ph<-c(0.5,0.3,0.2)
pxw<-matrix(c(0.4,0.3,0.2,0.1,0.15,0.1,0.35,0.4,0.2,0.4,0.15,0.25),nrow=H0,ncol=L,byrow=T)
  
K0=30
nk0<-16
n0=K0*nk0

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

components0<-rep(0,K0)
while(length(which(components0==1))!=round(K0*ph[1]) || length(which(components0==2))!=round(K0*ph[2]) || length(which(components0==3))!=round(K0*ph[3])){
  components0 <- sample(1:H0,prob=ph,size=K0,replace=TRUE)      
}

dati <- matrix(nrow=n0,ncol=2) 
dati[,1] <- seq(1:n0)
d1<-NULL
for(k in 1:K0){
  d1<-c(d1,rep(k,nk0))
}

dati[,2]<-d1

colh<-NULL
for(i in 1:(K0)){
  colh<-c(colh,rep(components0[i],nk0))
}
datic<-cbind(dati,colh)

daticc<-datic[order(datic[,3]),]
dati<-daticc[,1:2]

conteggio0<-c(length(which(components0==1))*nk0,length(which(components0==2))*nk0,length(which(components0==3))*nk0)

dati2<-NULL
samples2 <- NULL
for(j in (1:H0)){
  samples <- sample(1:L,prob=pxw[j,],size=conteggio0[j],replace=TRUE) 
  for(i in (1:conteggio0[j])){
   dati2p = cbind(sample(0:(c1-1),prob=p1[,samples[i]],size=1),sample(0:(c1-1),prob=p2[,samples[i]],size=1),sample(0:(c2-1),prob=p3[,samples[i]],size=1),sample(0:(c2-1),prob=p4[,samples[i]],size=1),sample(0:(c3-1),prob=p5[,samples[i]],size=1),sample(0:(c3-1),prob=p6[,samples[i]],size=1))
   dati2=rbind(dati2,dati2p)
  } 
  samples2<-c(samples2,samples)  
}

dati<-cbind(dati,dati2)  
datic<-cbind(dati,samples2)

s<-data
s<-as.matrix(s)

y<-list(s[1,])
for(i in (2:dim(s)[1])){
  y[[i]]<-s[i,]
}

main2(y)


}
