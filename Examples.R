setwd("./HawkesProcesses/paper/CodeForTimings")
library(Runuran)
library(microbenchmark)
source("HawkesProcessGeneration.R")


res<-PPsumExp.naive(lam=1.e4,a=3,b=4.5)
hist(res,breaks=60)
microbenchmark(PPsumExp.naive(lam=1.e4),times=100)


res<-PPsumExp(lam=1.e4,a=3,b=4.5)
hist(res,breaks=60)
microbenchmark(PPsumExp(lam=1.e4),times=100)


res<-PP.Algo1(lam=1.e4,a=3,b=4.5)
hist(res,breaks=60)
microbenchmark(PP.Algo1(lam=1.e5),times=100)


res<-PP.Algo2(lam=1.e4,a=3,b=4.5)
hist(res,breaks=60)
microbenchmark(PP.Algo2(lam=1.e5))# unsorted
microbenchmark(PP.Algo2(lam=1.e5,sortYN=TRUE))# sorted


res<- NHPP.Algo3_2(invCDF=qexp,Aphi=1.e4,sortYN=FALSE)
hist(res,breaks=60)
microbenchmark(NHPP.Algo3_2(invCDF=qexp,Aphi=1.e5,sortYN=FALSE))


res<- NHPP.Algo3_1(invCDF=qexp,Aphi=1.e4)
hist(res,breaks=60)
microbenchmark(NHPP.Algo3_1(invCDF=qexp,Aphi=1.e5))



res<- NHPPthinning.naive(phi=function(x) exp(-x)*1.e3, a=0,b=10)
hist(res$events,breaks=60,main=paste("Histogram, ",round(res$percAcc*100,2),"% were accepted"))
microbenchmark(NHPPthinning.naive(phi=function(x) exp(-x)*1.e4, a=0,b=10),times=20)


 res<- NHPP.Algo4(phi=function(x) exp(-x)*1.e4,UBnTr=11000, a=0.1,b=7)
hist(res$events,breaks=60,main=paste("Histogram Algo4, ",round(res$percAcc*100,2),"% were accepted"))
c(min(res$events),max(res$events))# 0.1000096, 6.930091
microbenchmark(NHPP.Algo4(phi=function(x) exp(-x)*1.e4, UBnTr=11000, a=0,b=10),times=20)
 res<-HP.Algo5(lam0=function(x)ifelse(x<1,100,0), phi=function(x) exp(-x)*0.5, a=0,b=100)
 c(length(res$events),res$nTr,res$percAcc) #[1] 199.000000 210.000000   0.947619
 hist(res$events,breaks=100,main=paste("Histogram Algo4, ",round(res$percAcc*100,2),"% were accepted"))
system.time(res<-HP.Algo5(lam0=function(x)ifelse(x<1,5.e3,0), phi=function(x) exp(-x)*0.5, a=0,b=100))
hist(res$events,breaks=100,main=paste("Histogram Algo4, ",round(res$percAcc*100,2),"% were accepted"))

res<- HP.Algo6(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=FALSE)
microbenchmark(res<- HP.Algo6(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=FALSE))
microbenchmark(res<- HP.Algo6(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=TRUE))
 hist(res$events,breaks=100)
 c(length(res$events),max(res$events),res$nRej,res$kGener)

res<- HP.Algo6_7(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=FALSE)
 hist(res$events,breaks=100)
c(length(res$events),max(res$events),res$nRej,res$kGener) 
#  [1] 99532.00000    14.02002       0.00000   8.00000
microbenchmark(res<- HP.Algo6_7(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=FALSE))
microbenchmark(res<- HP.Algo6_7(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=TRUE))


res<-PP.Algo7(m=1.e5,T=1.5,lam=1)
hist(res$events,breaks=60)
microbenchmark(PP.Algo7(m=1.e5,T=1.5,lam=1))

res<-PP.Algo7(m=1.e5,T=1.5,lam=1,sortYN=TRUE)
resInd<-getIndices(res)
F <-res$events[resInd$endInd] # final event of all simulated non-empty realisations
hist(F,breaks=60)

res<- NHPP.Algo8(m=2.5e4,invCDF=qexp,Aphi=4,sortYN=FALSE)
hist(res$events,breaks=100)
microbenchmark(res<- NHPP.Algo8(m=2.5e4,invCDF=qexp,Aphi=4,sortYN=FALSE))


res <- HP.Algo9(m=1000,lam0InvCDF=function(x)10*x,Alam0=50,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=FALSE)
hist(res$events,breaks=100)
length(res$events)#[1] 99312
microbenchmark(res <- HP.Algo9(m=1000,lam0InvCDF=function(x)10*x,Alam0=50,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=FALSE))

 res <- HP.Algo9(m=1.e5,lam0InvCDF=NULL,Alam0=0,immigrant=0,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=FALSE)
 hist(res$events,breaks=100)
 length(res$events)#[1] 100155 
 max(res$events) #[1] 23.96439
min(res$events) #[1] [1] 1.549722e-06  # shows that the 0 immigrants are not returned

res.immi <- rRenewalProc(m=1.e4,T=10.5,phiInvCDF=function(u) 0.9+0.2*u,sortYN=TRUE)#inter event times are U(0.9,1.1)
res <- HP.Algo9(m=1.e5,lam0InvCDF=NULL,Alam0=0,immigrant=res.immi,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=FALSE)
c(length(res$events),min(res$events))#[1] 1.990700e+05 9.000066e-01
hist(res$events,breaks=200,main="Renewal Hawkes Process")
microbenchmark( { res.immi <- rRenewalProc(m=5.e3,T=10.5,phiInvCDF=function(u) 0.9+0.2*u,sortYN=TRUE)# inter event times are U(0.9,1.1);
res <- HP.Algo9(m=1.e5,lam0InvCDF=NULL,Alam0=0,immigrant=res.immi,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=FALSE) } )


res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,markInvCDF=function(u) qgamma(u,shape=4,rate=4))
microbenchmark(res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,markInvCDF=function(u) qgamma(u,shape=4,rate=4)),times=20)
length(res$events)

rNo <- function(n)rnbinom(n,size=4,prob=1/(1+0.5/4)) # to get n*=0.5 we use:  prob = 1/(1+n*/size)
res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,rNoffspring=rNo,sortYN=TRUE)
length(res$events)
microbenchmark(res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),repr=0.5,rNoffspring=rNo,sortYN=TRUE),times=20)

oRun.offspr<-Runuran::dgt.new(dnbinom(0:50,size=4,prob=1/(1+0.5/4)),from=0)
rNoFast <- function(n){ Runuran::ur(oRun.offspr,n) }
microbenchmark(res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,rNoffspring=rNoFast,sortYN=TRUE))
length(res$events)#[1] 99692

res.immi <- PP.Algo7(m=50,T=10,lam=100,sortYN=TRUE) 
res <- markHPcluster(m=1.e5,T=Inf,immigrant=res.immi,phiInvCDF=qexp,repr=0.5,markInvCDF=qexp)
microbenchmark({res.immi <- PP.Algo7(m=50,T=10,lam=100,sortYN=TRUE) 
res <- markHPcluster(m=1.e5,T=Inf,immigrant=res.immi,phiInvCDF=qexp,repr=0.5,markInvCDF=qexp)},times =100)

oRun.offspr<-Runuran::dgt.new(dnbinom(0:50,size=1,prob=1/(1+0.5/1)),from=0)
rNoFast <- function(n){ Runuran::ur(oRun.offspr,n) }
res.immi <- PP.Algo7(m=50,T=10,lam=100,sortYN=TRUE) 
res <- markHPcluster(m=1.e5,T=Inf,immigrant=res.immi,phiInvCDF=qexp,rNoffspring=rNoFast)
microbenchmark( { res.immi <- PP.Algo7(m=50,T=10,lam=100,sortYN=TRUE) 
 res <- markHPcluster(m=1.e5,T=Inf,immigrant=res.immi,phiInvCDF=qexp,rNoffspring=rNoFast) } )


