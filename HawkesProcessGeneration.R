#############################################################################
PPsumExp.naive <- function(a=0,b=1,lam=1){
# generates sorted PP : on (a,b)  using paper Algo 1 "Sum of Exponentials"
# (a,b) ... ibterval in which the events are generated
# lam ... rate
 En <- (b-a)*lam 
 tv<- numeric(En+5*sqrt(En))
 t<-a; i<-1  # 1:
 repeat({   # 2:
   E <- rexp(1)# 3: Generate a standard exponential E.
   t <- t+E/lam; tv[i]<-t # 4: Set t <- t+E/lam and ti<-t.
   if(t > b){ # 6:
      return(tv[1:(i-1)]) }  # 7: return (t1 t2, ... t[i-1])
   i <- i+1 # 8:
}) # # 9:
}# end PPsumExp.naive()

#res<-PPsumExp.naive(lam=1.e4,a=3,b=4.5)
#hist(res,breaks=60)
#microbenchmark(PPsumExp.naive(lam=1.e4),times=100)
#                        expr    min     lq     mean median     uq    max neval
# PPsumExp.naive(lam = 10000) 6.4449 6.8192 7.065276 6.9671 7.1828 8.3252   100

#############################################################################
PPsumExp <- function(a=0,b=1,lam=1){
# generates sorted PP on (a,b)  using paper Algo 1 "Sum of Exponentials"
# lam ... rate
# Note: uses a vector to pregenerate the exponential RVs; thus about 4x faster than PPsumExp.naive
 En <- (b-a)*lam
 tv<- numeric(En+5*sqrt(En)); Ev<- rexp(En+5*sqrt(En))
 t<-a; i<-1  #  Set t <- a and i <- 1.
 repeat({
   E <- Ev[i]#Generate a standard exponential E.
   t <- t+E/lam; tv[i]<-t # Set t <- t+E/(lam (b- a))and ti<-t.
#   print(paste("t=",t,"b=",b))
   if(t > b){ return(tv[1:(i-1)]) }   #return (t1 t2, ... t[i-1])
   i <- i+1
})
}# end PPsumExp()

#res<-PPsumExp(lam=1.e4,a=3,b=4.5)
#hist(res,breaks=60)
#microbenchmark(PPsumExp(lam=1.e4),times=100)
#                  expr    min      lq     mean  median      uq    max neval
# PPsumExp(lam = 10000) 1.4473 1.49505 1.557389 1.51305 1.55075 2.0566   100

#############################################################################
PP.Algo1 <- function(a=0,b=1,lam=1){
# uses cumsum() and a simple upper bound for the maximal length, to speed up; 
# generates sorted events of a Poisson Process on (a,b)  using 
# Algo 1 "Sum of Exponentials"
# lam ... rate
 En <- (b-a)*lam 
 tv<- cumsum(c(a,rexp(En+5*sqrt(En))/lam))[-1]
 return(tv[tv < b])
}# end PP.Algo1()

#res<-PP.Algo1(lam=1.e4,a=3,b=4.5)
#hist(res,breaks=60)
#microbenchmark(PP.Algo1(lam=1.e5),times=100)
#                         expr   min     lq     mean median     uq     max neval
# PP.Algo1(lam = 1e+05) 4.1754 4.862 5.109356  4.946 5.0222 15.657  1000


#############################################################################
PP.Algo2 <- function(a=0,b=1,lam=1,sortYN=FALSE){
#PPunifOrdStat<- function(a=0,b=1,lam=1,sortYN=FALSE){
# generates unsorted events of a Poisson Process on (a,b)  using
# Algo 2 " uniform order statistics"
# (a,b) ... ibterval in which the events are generated
# lam ... rate
# sortYN ... use TRUE if a sorted output is required
 N <- rpois(1,lam*(b-a)) # 1:
 U <- a+(b-a)*runif(N)   # 2:
 if(!sortYN) return(U)   
 return(sort(U))         # 3: 
}# end PP.Algo2()

#res<-PP.Algo2(lam=1.e4,a=3,b=4.5)
#hist(res,breaks=60)
#microbenchmark(PP.Algo2(lam=1.e5))# unsorted
#                  expr    min     lq     mean median     uq    max neval
# PP.Algo2(lam = 1e+05) 1.5018 1.5618 1.649186 1.5843 1.6082 5.5296   100
#microbenchmark(PP.Algo2(lam=1.e5,sortYN=TRUE))# sorted
#                                 expr    min     lq     mean  median     uq
# PP.Algo2(lam = 1e+05, sortYN = TRUE) 6.0066 6.3055 6.537556 6.36705 6.4971

####################################################################################
NHPP.Algo3_2 <- function(invCDF=qexp,Aphi=1.e4,sortYN=FALSE){
# generates events of NHPP using paper Algo 3 "NHPP Inversion" and 
#      Algo 2 (orderstatistic) to generate the Poisson Process
# invCDF ... inverse CDF of the intensity function normed to integral 1
# Aphi ... area below the intensity function; equal to the expected number of events
# sortYN ... use TRUE if a sorted output is required
 N <- rpois(1,Aphi)  # step 2.1
 V <- runif(N)   # step 2.2
 if(!sortYN) return(invCDF(V)) # 3:
 return(invCDF(sort(V)))
}# end NHPP.Algo3_2()
#res<- NHPP.Algo3_2(invCDF=qexp,Aphi=1.e4,sortYN=FALSE)
#hist(res,breaks=60)
#microbenchmark(NHPP.Algo3_2(invCDF=qexp,Aphi=1.e5,sortYN=FALSE))
#    min      lq     mean  median     uq    max neval
# 4.4245 4.77065 4.931673 4.84555 4.9554 8.5035   100

####################################################################################
NHPP.Algo3_1 <- function(invCDF=qexp,Aphi=1.e4){
# generates sorted events of a NHPP using paper Algo 3 "NHPP Inversion" and 
#      Algo 1 (sum exponentials) to generate the Poisson Process
# invCDF ... inverse CDF of the intensity function normed to integral 1
# Aphi ... area below the intensity function; equal to the expected number of events
 tv<- cumsum(c(0,rexp(Aphi+5*sqrt(Aphi))/Aphi))[-1]
 return(invCDF(tv[tv < 1]))
}# end NHPP.Algo3_1()
#res<- NHPP.Algo3_1(invCDF=qexp,Aphi=1.e4)
#hist(res,breaks=60)
#microbenchmark(NHPP.Algo3_1(invCDF=qexp,Aphi=1.e5))
#     lq     mean median      uq     max neval
# 7.8278 8.621239 8.4391 8.73115 12.4151   100


####################################################################################
NHPPthinning.naive<- function(phi=function(x) exp(-x)*1.e3, a=0,b=10){
# generates sorted events of NHPP on (a,b) using a slow implementation of
# Algo 4 "Thinning with modified upper bound" for decreasing intensity phi
# returns list with the vector events and percAcc= proportion of accepted events
# phi ... intensity function, area below is expected number of events

#phi=function(x) exp(-x); a=0;b=10
 t<-a; j<-0; phit<- phi(a); tv <- NULL  # 1
 repeat({                #2 
   t <- t + rexp(1)/phit # 3
   if(t > b){return(list(events=tv,percAcc=length(tv)/j))}# 4 & 5
   U<- runif(1)*phit; phit <- phi(t); # 7
   if(U <= phit){        #8
      tv <- c(tv,t);  # 9
   }
   j <- j+1
})
}   
#res<- NHPPthinning.naive(phi=function(x) exp(-x)*1.e3, a=0,b=10)
#hist(res$events,breaks=60,main=paste("Histogram, ",round(res$percAcc*100,2),"% were accepted"))
#microbenchmark(NHPPthinning.naive(phi=function(x) exp(-x)*1.e4, a=0,b=10),times=20)
#      min       lq     mean   median       uq      max neval
# 102.2551 107.4632 108.9492 109.2177 110.2408 120.6769    20


##########################################
NHPP.Algo4 <- function(phi=function(x) exp(-x)*1.e3,UBnTr, a=0,b=10){
# generates sorted events of NHPP on (a,b) using a faster implementation of
# Algo 4 "Thinning with modified upper bound" for decreasing intensity phi
# returns list with   events = event-vector, nTr= number of Trials, 
#                     and percAcc= percentaage of accepted events
# phi ... intensity function
# UBnTr... maximal number of events tried; 
#          for generating all events UBnTr must be larger than area below phi
# (a,b) ... ibterval in which the events are generated
 t<-a; i<-1; phit<- phi(a); tv <- numeric(UBnTr)  # 1
 ruv <- runif(UBnTr); rexpv <- rexp(UBnTr)
 for(j in 1:UBnTr){
   t <- t + rexpv[j]/phit # 3
   if(t > b){ if(j==1){ return(list(events=NULL,nTR=0,percAcc=NA)) } # 4 & 5
      return(list(events=tv[1:(i-1)],nTr=(j-1),percAcc=(i-1)/(j-1)))}
   U<- ruv[j]*phit; phit <- phi(t); # 7
   if(U <= phit){         # 8
      tv[i] <- t; i<-i+1  # 9
   }
 }
}
# res<- NHPP.Algo4(phi=function(x) exp(-x)*1.e4,UBnTr=11000, a=0.1,b=7)
#hist(res$events,breaks=60,main=paste("Histogram Algo4, ",round(res$percAcc*100,2),"% were accepted"))
#c(min(res$events),max(res$events))# 0.1000096, 6.930091
#microbenchmark(NHPP.Algo4(phi=function(x) exp(-x)*1.e4, UBnTr=11000, a=0,b=10),times=20)
#    min     lq    mean  median     uq    max neval
# 5.8393 6.0711 6.45611 6.24355 6.6933 7.8268    20


###########################################################
HP.Algo5 <- function(lam0=function(x)10, a=0,b=10, phi=function(x) exp(-x)*0.5,maxTrial=50000){
# generates sorted events of a Hawkes Process on (a,b)  using paper 
# Algo 5 "Ogata's Modified Thinning" for decreasing phi
# returns list with   events = event-vector, nTr= number of Trials, 
#                     and percAcc= percentaage of accepted events
# lam0 ... base intensity function on (a,b)
# (a,b) ... ibterval in which the events are generated
# phi ... excitation function
# maxTrial ... maximal number of Trials

 BnTr=100 # number of random variates that are precalcuated; 100 should not be changed
 t<-a; lamt<-lam0(a); phi0<- phi(0); i<-1; # 1:
 tv <- NULL 
 for(j in 1:maxTrial){ # 2:
   if(j%%BnTr==1){
    ruv<- runif(BnTr); rexpv <- rexp(BnTr); tv <- c(tv,numeric(BnTr))# pre calculation of random variates
   }
   t <- t + rexpv[1+j%%BnTr]/lamt          # 3:
   if(t > b){                              # 4:
      if(j==1){ return(list(events=NULL,nTR=0,percAcc=NA)) }
      return(list(events=tv[1:(i-1)],nTr=(j-1),percAcc=(i-1)/(j-1)) )} # # 5:
   U<- ruv[1+j%%BnTr]*lamt;                # 7: 1 
   lamt <- lam0(t)+sum(phi(t-tv[1:j]));    # 7: 2
   if(U <= lamt){                          # 8:
      tv[i] <- t; i<-i+1; lamt<-lamt+phi0  # 9:
   }# 10:
 }# 11:
  return(list(events=tv[1:(i-1)],nTr=(j)))
}# end HP.Algo5()
# res<-HP.Algo5(lam0=function(x)ifelse(x<1,100,0), phi=function(x) exp(-x)*0.5, a=0,b=100)
# c(length(res$events),res$nTr,res$percAcc) #[1] 199.000000 210.000000   0.947619
# hist(res$events,breaks=100,main=paste("Histogram Algo4, ",round(res$percAcc*100,2),"% were accepted"))
#system.time(res<-HP.Algo5(lam0=function(x)ifelse(x<1,5.e3,0), phi=function(x) exp(-x)*0.5, a=0,b=100))
#   user  system elapsed 
#   2.21    0.00    2.20 
#hist(res$events,breaks=100,main=paste("Histogram Algo4, ",round(res$percAcc*100,2),"% were accepted"))



################################################################################
#################################################################################

HP.Algo6 <- function(lam0InvCDF=function(x)x,Alam0=5,phiInvCDF=qexp,repr=0.5,T=100,sortYN=FALSE){
# generates events of a Hawkes Process on (0,T)  using paper Algo 6 "Cluster based".
# returns list:   events = event-vector, nRej= # of rejected event, kGener= # of generations
#                 note that events are only rejected if they are > T  !
# lam0invCDF ... inverse CDF of the base intensity function normed to integral 1
#                it is used to generate the immigrants
# Alam0 ... area below the base intensity function; equal to the expected number of immigrants
# sortYN ... use TRUE if a sorted output is required
# phiInvCDF ... inverse CDF of the excitation function normed to integral 1
# repr=0.5 ... reproduction number n* = expected number of offsprings
# T ... Inf or finite ... upper bound for the simulation
# sortYN ... use TRUE if a sorted output is required
 G <- list()
 N <- rpois(1,Alam0)  # 1.1:
 if(N==0) return(NULL)
 nRej <- 0
 G[[1]] <- lam0InvCDF(runif(N))# 1.2:
 for(k in 2:100000){           # 2:
   N <- rpois(length(G[[k-1]]),repr)  # 3:&4:
   if(sum(N)>0){                      # 5:
     G[[k]] <- rep(G[[k-1]],N) + phiInvCDF(runif(sum(N)))# 6:&7:&8:
   }else{                                                # 12.1:
      if(!sortYN){return(list(events=unlist(G),nRej=nRej,kGener=k))} 
	  return(list(events=sort(unlist(G)),nRej=nRej,kGener=k)) # 12.2
   }	  
   if(T<Inf){ 
     iv <-  ( G[[k]] < T) # 9:
	 nRej <- nRej + sum(!iv)
     G[[k]] <- G[[k]][iv] # 9:
   }
 }
}# end HP.Algo6()  
#res<- HP.Algo6(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=FALSE)
#microbenchmark(res<- HP.Algo6(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=FALSE))
#    min      lq     mean median      uq   max neval
# 4.3423 4.40465 4.481202 4.4425 4.49325 5.477   100
#microbenchmark(res<- HP.Algo6(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=TRUE))
#    min      lq    mean  median      uq     max neval
# 9.2318 9.43985 10.0981 10.2071 10.4142 19.3562   100
# hist(res$events,breaks=100)
# c(res$events,max(res$events),res$nRej,res$kGener)



###############################################

HP.Algo6_7 <- function(lam0InvCDF=function(x)x,Alam0=5,phiInvCDF=qexp,repr=0.5,T=100,sortYN=FALSE){
# generates events of a Hawkes Process on (0,T)  using paper Algo 6 "Cluster based" together
#         with Algo 7 to generate the NHPP-subprocesses of the Hawkes process.
# returns list:   events = event-vector, nRej= # of rejected event, kGener= # of generations
#                 note that events are only rejected if they are > T  !
# lam0invCDF ... inverse CDF of the base intensity function normed to integral 1
#                it is used to generate the immigrants
# Alam0 ... area below the base intensity function; equal to the expected number of immigrants
# sortYN ... use TRUE if a sorted output is required
# phiInvCDF ... inverse CDF of the excitation function normed to integral 1
# repr=0.5 ... reproduction number n* = expected number of offsprings
# T ... Inf or finite ... upper bound for the simulation
# sortYN ... use TRUE if a sorted output is required
 G <- list()
 N <- rpois(1,Alam0)  # step 1.1:
 if(N==0) return(NULL)
 nRej <- 0
 G[[1]] <- lam0InvCDF(runif(N))   # step 1.2:
 for(k in 2:100000){
   m <- length(G[[k-1]])
   N <- rpois(1,m*repr)  # Algo 7 1:
   if(N>0){	    
     Um <- runif(N)*m    # Algo 7 2:
     IDnoM1 <- as.integer(Um)  # Algo 7 3:
     G[[k]] <- G[[k-1]][1+IDnoM1] + phiInvCDF(Um-IDnoM1) # Algo 6 6:&7:&8:   and Algo 7 4:
   }else{
      if(!sortYN){return(list(events=unlist(G),nRej=nRej,kGener=k))}
	  return(return(list(events=sort(unlist(G)),nRej=nRej,kGener=k)))
   }	  
   if(T<Inf){ 
     iv <-  ( G[[k]] < T)
	 nRej <- nRej + sum(!iv)
     G[[k]] <- G[[k]][iv] 
   }
 }
}# end HP.Algo6_7()   
#res<- HP.Algo6_7(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=FALSE)
# hist(res$events,breaks=100)
#c(length(res$events),max(res$events),res$nRej,res$kGener) 
#  [1] 99532.00000    14.02002       0.00000   8.00000
#microbenchmark(res<- HP.Algo6_7(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=FALSE))
#    min      lq     mean  median     uq    max neval
# 2.4928 2.55075 2.621869 2.57185 2.5992 3.5936   100
#microbenchmark(res<- HP.Algo6_7(lam0InvCDF=function(x)x,Alam0=80000,phiInvCDF=qexp,repr=0.2,T=Inf,sortYN=TRUE))
#    min      lq    mean  median     uq     max neval
# 7.4712 8.0612 8.370767 8.25255 8.3365 24.3811   100



#########################################################################################################################
#########################################################################################################################
#   start of code to generate m independent realisations
#########################################################################################################################
getIndices <- function(res){
# extrats from the list returned by the algorithms below with sortYN =TRUE 
# for each simulated independent realisation the begin and end location in the events vector.
# res ... list returned by rHawkesMarks() or rHawkes() .
# TODO use the unsorted() function to check if the list is sorted, otherwise sort !!!!!!!!!!!!!!!!!!
  if(is.unsorted(res$IDno)){
    print("getIndices() works only for a list with sorted events.!!!")
    return(NULL)	
  }
  wh1 <- which(0!=diff(res$IDno))
  endInd <- c(wh1,length(res$IDno)) # allows eg. to get easily a vector of the largest event of all repetitions
  startInd <- c(1,wh1+1)
  lengths <- endInd-startInd +1
  return(list(startInd=startInd,endInd=endInd))
}

###############################################

PP.Algo7 <- function(m,T=1,lam=1,sortYN=FALSE){
# generates m inependent realisation of a Poisson Process on (0,T)  
# using paper Algo 7 "uniform order statistics"
# returns a list holding events = event-vector, IDno = vector of the ID numbers
# m ... number of independent repetitions
# T ... upper border of domain
# lam ... rate
# sortYN ... use TRUE if a sorted output is required
 N <- rpois(1,T*lam*m)    # 1:
 U <- runif(N)*m       # 2:
 IDnoM1 <- as.integer(U) # 3.1:
 if(!sortYN) return( list(events=T*(U-IDnoM1),IDno=IDnoM1+1) ) # 3.2 & 4:
 IDno <- IDnoM1+1
 U <- U-IDnoM1
 oIDev<-order(IDno,U)
 return(list(events=T*U[oIDev],IDno=IDno[oIDev]))
}# end PP.Algo7()   
#res<-PP.Algo7(m=1.e5,T=1.5,lam=1)
#hist(res$events,breaks=60)
#microbenchmark(PP.Algo7(m=1.e5,T=1.5,lam=1))
#   min      lq     mean median     uq    max neval
#3.0722 3.18235 3.329134  3.264 3.3495 4.5935   100

#res<-PP.Algo7(m=1.e5,T=1.5,lam=1,sortYN=TRUE)
#resInd<-getIndices(res)
#F <-res$events[resInd$endInd] # final event of all simulated non-empty realisations
#hist(F,breaks=60)

##############################################################
NHPP.Algo8 <- function(m=1,invCDF=qexp,Aphi=2,sortYN=FALSE){
# simuates m independent repetitons of NHPP using paper Algo 3 "NHPP Inversion" 
# and Algo 7 to generate the Poisson processes. 
# returns a list holding events = event-vector, IDno = vector of the ID numbers
# m ... number of independent repetitions that should be generated
# invCDF ... inverse CDF of the intensity function normed to integral 1
# Aphi ... area below the intensity function; equal to the expected number of events
# sortYN ... use TRUE if a sorted output is required

 PP <- PP.Algo7(lam=Aphi,m=m,sortYN=sortYN) # 1
 PP$events <- invCDF(PP$events)             # 2.1
 return(PP)                                 # 2.2
}# end PP.Algo8()   
#res<- NHPP.Algo8(m=2.5e4,invCDF=qexp,Aphi=4,sortYN=FALSE)
#hist(res$events,breaks=100)
#microbenchmark(res<- NHPP.Algo8(m=2.5e4,invCDF=qexp,Aphi=4,sortYN=FALSE))
#   min      lq     mean  median      uq   max neval
# 5.064 5.17975 5.636845 5.24385 6.11545 8.033   100

##########################################
##########################################
##########################################
##########################################
HP.Algo9 <- function(m=1,lam0InvCDF=NULL,Alam0=5,immigrant=0,phiInvCDF=qexp,repr=0.5,T=100,sortYN=FALSE){
# simulates m independent realisations of HP on (0,T) using Algorithm 9 (cluster based method) 
#         together with Algo 7 to generate the NHPP-subprocesses of the Hawkes process.
# returns list:   events = event-vector, IDno = IDno vector,
#                   nRej = # of rejected event, kGener = # of generations
#                 note that events are only rejected if they are > T  !
# m ... number of independent realisations that are generated
#       m is ignored if lam0InvCDF is not given and immigrant is list of simulated immigrant events
# lam0InvCDF ... inverse CDF of the base intensity function normed to integral 1
#                it is used to generate the immigrants
# Alam0 ... area below the base intensity function; equal to the expected number of immigrants
#           Alam0 is ignored if lam0InvCDF is not a function.
# immigrant ... vector of fixed immigrant events on (0,T); is used for all realisations and thus repeated n times
#                    the values of these "fixed immigrants" are NOT returned in the final list that is returned 
#               or list holding m simulated realisations of an immigrant process on (0,T).
#                      these simulated immigrant values are also included in the finally returned list            
#           immigrant is ignored if lam0InvCDF is a function.
# phiInvCDF ... inverse CDF of the excitation function normed to integral 1
# repr=0.5 ... reproduction number n* = expected number of offsprings
# T ... Inf or finite ... upper bound for the simulation
# sortYN ... use TRUE if a sorted output is required

memBound=2^28 # to stop a simulation that requires too much memory
              # necessary esepcially for experiments with repr very close or equal or larger than 1  !!!!!

 nRej <- 0
 G <- ID <- list()
 if(is.function(lam0InvCDF)){
   immigrant <- NHPP.Algo8(m=m,invCDF=lam0InvCDF,Aphi=Alam0,sortYN=FALSE) # 1:
   G[[1]] <- immigrant$events
   ID[[1]] <- immigrant$IDno
 }else if(is.list(immigrant)){
   G[[1]] <- immigrant$events
   ID[[1]] <- immigrant$IDno
 }else{
   G[[1]] <- rep(immigrant,m)
   ID[[1]] <- rep(1:m,each=length(immigrant))
 }
#print(rbind(round(G[[1]],3),ID[[1]]))
 for(k in 2:100000){                                                    # 2:
#print(paste("k=",k))
   lkm1 <- length(G[[k-1]])
#print(lkm1)
#print(sortYN) 
   offspr <- NHPP.Algo8(m=lkm1,invCDF=phiInvCDF,Aphi=repr,sortYN=FALSE) # 3
#print(rbind(round(offspr$events,3),offspr$IDno))
   N <- length(offspr$events)
   if(N>0){	                                              # 4: condition negated
    if(sum(N)*8+object.size(G)>memBound) {
	    print("total memor size required larger than memBound = 250 MB")
	    print("error in rHawkes; object too large; reduce repr,T or n!!")
        return(NULL)
   	 }	

     ID[[k]] <- ID[[k-1]][offspr$IDno]                 # 7
     G[[k]] <- G[[k-1]][offspr$IDno] + offspr$events   # 8
#print(rbind(round(G[[k]],3),ID[[k]]))
   }else{
      if( !is.function(lam0InvCDF) & ! is.list(immigrant)){# "fixed immigrants" not included in the final list
	      G <- G[-1]; ID<-ID[-1] }
      if(!sortYN){return(list(events=unlist(G),IDno=unlist(ID),nRej=nRej,kGener=k-1))} # 5
	  G <- unlist(G); ID <- unlist(ID)
	  oIDev<-order(ID,G)
      return(list(events=G[oIDev],IDno=ID[oIDev],nRej=nRej,kGener=k-1))	  
   }                         # 10:	  
   if(T<Inf){ 
     iv <-  ( G[[k]] < T)    # 9:
	 nRej <- nRej + sum(!iv) # 9:
     G[[k]] <- G[[k]][iv]    # 9:
     ID[[k]] <- ID[[k]][iv]  # 9:
   }
 }       # 11:
}# end HP.Algo9()     

# example for immigrants simulated as NHPP uing "lam0InvCDF" and "Alam0"
#res <- HP.Algo9(m=1000,lam0InvCDF=function(x)10*x,Alam0=50,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=F)
#hist(res$events,breaks=100)
#length(res$events)#[1] 99312
#microbenchmark(res <- HP.Algo9(m=1000,lam0InvCDF=function(x)10*x,Alam0=50,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=F))
#    min     lq     mean median     uq    max neval
# 4.9244 5.2123 5.387533 5.3261 5.4118 7.4739   100

# example: simulate a "Cluster" a HP with a single immigrant at 0 (it is not stored in the result!!) 
# res <- HP.Algo9(m=1.e5,lam0InvCDF=NULL,Alam0=0,immigrant=0,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=FALSE)
# hist(res$events,breaks=100)
# length(res$events)#[1] 100155 
# max(res$events) #[1] 23.96439
# min(res$events) #[1] [1] 1.549722e-06  # shows that the 0 immigrants are not returned

# example: Simulate HP with Renewal-Process immigrants (called Renewal Hawkes Process)
# first we generate the immigrants
# res.immi <- rRenewalProc(m=1.e4,T=10.5,phiInvCDF=function(u) 0.9+0.2*u,sortYN=TRUE)# inter event times are U(0.9,1.1)
# res <- HP.Algo9(m=1.e5,lam0InvCDF=NULL,Alam0=0,immigrant=res.immi,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=FALSE)
# c(length(res$events),min(res$events))#[1] 1.990700e+05 9.000066e-01
# hist(res$events,breaks=200,main="Renewal Hawkes Process")
# microbenchmark( { res.immi <- rRenewalProc(m=5.e3,T=10.5,phiInvCDF=function(u) 0.9+0.2*u,sortYN=TRUE)# inter event times are U(0.9,1.1);
# res <- HP.Algo9(m=1.e5,lam0InvCDF=NULL,Alam0=0,immigrant=res.immi,phiInvCDF=qexp,repr=0.5,T=Inf,sortYN=FALSE) } )
#   min     lq     mean  median     uq    max neval
# 5.602 5.8801 6.067201 5.98705 6.1096 7.2466   100



####################################################################################################
####################################################################################################
####################################################################################################
markHPcluster<-function(m=1,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,markInvCDF=qexp, 
                         rNoffspring=NULL,sortYN=TRUE){
# Generates events of n independent Hawkes processes with iid marks using the 
# cluster based method with rejection events larger T in every generation.
# returns a list with the first vector holding all events and the 
# second vector the ID-numbers of these events.
# 
# m ... number of indepednent realisations that should be generated; 
#       m is ignored if immigrant is list of simulated immigrant events
# T ... upper border of event times to generate
# immigrant ... vector of fixed immigrant events on (0,T); is used for all realisations and thus repeated n times
#                    the values of these "fixed immigrants" are NOT returned in the final list that is returned 
#               or list holding m simulated realisations of an immigrant process on (0,T).
#                      these simulated immigrant values are also included in the finally returned list            
# phiInvCDF ... inverse CDF of self exciting density phi; need not have integral 1 
# repr ... reproduction number (mean number of off-springs); 
#          repr is ignored if dNoffspringV is set !!
# markInvCDF ... inverse CDF of distribution of the marks Yi. When the marks-distribution
#            has expectation 1 repr can be used to select the average reproduction number.
#            markInvCDF is ignored when dNoffspringV is set !!
# rNoffspring ... generates RVs of the distribution of the number of offsprings which is 
#                  typically a Poisson mixture distribution.
# sortYN ... use TRUE if a sorted output is required
memBound=2^28

 nRej <- 0
 G <- ID <- list()
 if(is.list(immigrant)){
   G[[1]] <- immigrant$events
   ID[[1]] <- immigrant$IDno
 }else{
   G[[1]] <- rep(immigrant,m)
   ID[[1]] <- rep(1:m,each=length(immigrant))
 }
 for(k in 2:100000){ # 2:####################################
   lkm1 <- length(G[[k-1]])
   if(!is.function(rNoffspring)){
	  N <- rpois(lkm1,repr*markInvCDF(runif(n=lkm1))) # simulates marks to get expectation of Poisson
   }else{ 
       N <- rNoffspring(lkm1)
   }
   if(sum(N)>0){ # 5:
     if(sum(N)*8+object.size(G)>memBound) {
	    print("total memor size required larger than memBound = 250 MB")
	    print("error in rHawkes; object too large; reduce repr,T or n!!")
        return(NULL)
   	 }	
     G[[k]] <- rep(G[[k-1]],N) + phiInvCDF(runif(sum(N)))# 6:&7:&8:
     ID[[k]]  <- rep(ID[[k-1]], times=N)    
   }else{ 
     if(!is.list(immigrant)){G <- G[-1]; ID<-ID[-1] }# fixed immgrants are not returned in the final list
     if(!sortYN){return(list(events=unlist(G),IDno=unlist(ID),nRej=nRej,kGener=k-1))}
 	 G <- unlist(G); ID <- unlist(ID)
 	 oIDev<-order(ID,G)
     return(list(events=G[oIDev],IDno=ID[oIDev],nRej=nRej,kGener=k-1))	  
   }                         # 10:	  
   if(T<Inf){ 
     iv <-  ( G[[k]] < T) # 9:
	 nRej <- nRej + sum(!iv)
     G[[k]] <- G[[k]][iv] # 9:
   }
 }
}# end of markHPcluster()

########
# example: simulating marked HP with marks following Gamma shape=4 rate=4; E(Gamma)=1
# we use repr to select the value of n*=0.5
# with direct simulation of the marks we get (quite slow):
#res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,markInvCDF=function(u) qgamma(u,shape=4,rate=4))
#microbenchmark(res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,markInvCDF=function(u) qgamma(u,shape=4,rate=4)),times=20)
#length(res$events)
#      min       lq     mean   median       uq     max neval
# 177.9746 182.3557 186.2465 187.1231 187.9244 195.615    20

# some speed up is possible if we use the Poisson mxture distribution for gamma marks is negative binomial with size 4
#rNo <- function(n)rnbinom(n,size=4,prob=1/(1+0.5/4)) # to get n*=0.5 we use:  prob = 1/(1+n*/size)
#res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,rNoffspring=rNo,sortYN=TRUE)
#length(res$events)
#microbenchmark(res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),repr=0.5,rNoffspring=rNo,sortYN=TRUE),times=20)
#   min      lq     mean  median       uq     max neval
# 41.29 41.5586 44.81024 42.9154 46.56815 53.2822    20

# extra speedup if we use the guide table method to generate from the negative binomial distribution
#oRun.offspr<-Runuran::dgt.new(dnbinom(0:50,size=4,prob=1/(1+0.5/4)),from=0)
#rNoFast <- function(n){ Runuran::ur(oRun.offspr,n) }
#microbenchmark(res <- markHPcluster(m=1.e5,T=Inf,immigrant=c(0),phiInvCDF=qexp,repr=0.5,rNoffspring=rNoFast,sortYN=TRUE))
#     min       lq     mean  median       uq     max neval
# 15.2095 16.20695 17.41199 17.1409 18.00375 25.1997   100
#length(res$events)#[1] 99692

#####################
# another example of simulating a marked HP with immigrants following a PP on (0,10)
# and both excitation function and marks follow the exponential distribution.
#res.immi <- PP.Algo7(m=50,T=10,lam=100,sortYN=TRUE) 
#res <- markHPcluster(m=1.e5,T=Inf,immigrant=res.immi,phiInvCDF=qexp,repr=0.5,markInvCDF=qexp)
#microbenchmark({res.immi <- PP.Algo7(m=50,T=10,lam=100,sortYN=TRUE) 
#res <- markHPcluster(m=1.e5,T=Inf,immigrant=res.immi,phiInvCDF=qexp,repr=0.5,markInvCDF=qexp)},times =100)
#     min      lq     mean   median      uq     max neval
# 27.9828 28.6082 29.41321 28.84105 29.7344 36.7374   100

#oRun.offspr<-Runuran::dgt.new(dnbinom(0:50,size=1,prob=1/(1+0.5/1)),from=0)
#rNoFast <- function(n){ Runuran::ur(oRun.offspr,n) }
#res.immi <- PP.Algo7(m=50,T=10,lam=100,sortYN=TRUE) 
#res <- markHPcluster(m=1.e5,T=Inf,immigrant=res.immi,phiInvCDF=qexp,rNoffspring=rNoFast)
#microbenchmark( { res.immi <- PP.Algo7(m=50,T=10,lam=100,sortYN=TRUE) 
#  res <- markHPcluster(m=1.e5,T=Inf,immigrant=res.immi,phiInvCDF=qexp,rNoffspring=rNoFast) } )
#     min      lq     mean   median       uq     max neval
# 17.8354 19.5716 20.33421 19.83205 20.64185 33.2306   100


############################################################################################
############################################################################################
############################################################################################
rRenewalProc <- function(m=1,T,phiInvCDF=function(u) 0.9+0.2*u, sortYN=FALSE){
# Generates events for n independent Renewal Processes
# returns a list with the first vector holding all events and the 
# second vector the corresponding IDnumbers (ie. the number of the repetition.)
#
# m ... number of indepednent realisations that will be generated
# T ... upper border of event times to generate; T must be finite for renewal processes
# phiInvCDF ... inverse CDF of the inter-event times
# sortYN ... use TRUE if a sorted output is required

 evev <- phiInvCDF(runif(m))
 i<- 1; IDno<- events <- list()
 while(TRUE){
   iv <- which(evev<T)
   if(length(iv)==0) break()
   events[[i]] <- evev[iv]
   if(i==1){ IDno[[i]] <- iv
   }else{IDno[[i]] <- IDno[[i-1]][iv]}
   evev <- events[[i]] + phiInvCDF(runif(length(events[[i]])))
   i<- i+1
 }  
 if(!sortYN) return(list(events=unlist(events),IDno=unlist(IDno)   ))
 IDno=unlist(IDno)
 ordID <- order(IDno)# as the events are generated in increasing order
                     # it is enough to sort only with respect to IDno
 return(list(events=unlist(events)[ordID],IDno=IDno[ordID]))
}# end of rRenewalProc()

# res <- rRenewalProc(m=1.e4,T=10,phiInvCDF=function(u) 0.9+0.2*u)
# hist(res$events,breaks=100,main="renewal Process with inter event-times ~ U(0.9,1.1)")








