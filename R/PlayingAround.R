library(ape)
library(picante)
library(mvtnorm)
library(boot)
library(deSolve)

htree<-rtree(12)
ptree<-rtree(15)

hcor<-vcv(htree,corr=T)
pcor<-vcv(ptree,corr=T)

hdim<-dim(hcor)[1]
pdim<-dim(pcor)[1]

Jh<-matrix(rep(1,hdim^2),nrow=hdim,ncol=hdim)
Jp<-matrix(rep(1,pdim^2),nrow=pdim,ncol=pdim)

Ah<-hcor
Ap<-pcor

Ih<-diag(hdim)
Ip<-diag(pdim)

ele1<-kronecker(Jp,Ah) #hosts nested in parasites (first go through all hosts for parasite 1, then all hosts for parasite 2, etc.)
ele2<-kronecker(Ap,Jh)
ele3<-kronecker(Ip,Ah)
ele4<-kronecker(Ap,Ih)
ele5<-kronecker(Ap,Ah)

ele6<-kronecker(Jp,Ih)
ele7<-kronecker(Ip,Jh)
ele8<-kronecker(Ip,Ih)

a<-5
b1<-0.1
b2<-0.1
b3<-0.2
b4<-0.1
b5<-0.2
b6<-0.1
b7<-0.1
b8<-0.1

fullcor<-a*(b1*ele1+b2*ele2+b3*ele3+b4*ele4+b5*ele5+b6*ele6+b7*ele7+b8*ele8)

test<-rmvnorm(1,mean=rep(-1,180),sigma=fullcor)
tt<-rbinom(180,1,inv.logit(test))
tmat<-matrix(tt,nrow=pdim,ncol=hdim)

model<-function(t,y,params){
 #died<-which(rbinom(params$nsite,1,params$d)==1)
 siteocc<-sum(y)
 #print(y)
 if (siteocc>0){
  transprob <- y / siteocc
 } else {transprob <- rep(0,params$nspec)}
 inprob <- 1/params$nspec
 #print(siteocc)
 transprob[transprob<0]<-0
 #print(transprob)
 #diespec <- sample.int(params$nspec,siteocc,F,params$d*transprob)
 die<-rbinom(params$nspec,siteocc,params$d*transprob)
 #print(die)
 torep<-sum(die)+params$nsite-siteocc
 newy <- y - die
 if (torep>0){
  repspec <- sample.int(params$nspec,torep,T,(1-params$poolprob)*transprob+params$poolprob*inprob)
  tab <- table(repspec)
  newy[as.numeric(names(tab))] <- newy[as.numeric(names(tab))] + tab 
 } else {newy <- y}
 return(list(newy))
}


model2<-function(t,y,params) {
  empt<-y==0
  #print(empt)
  dead<-c(sample(which(!empt),rbinom(1,sum(!empt),params$d)),which(empt))
  #print(dead)
  specnums<-table(y)
  #print(specnums)
  specprobs2<-rep(0,params$nspec)
  nam<-as.numeric(names(specnums))
  if (0 %in% nam){nums<-specnums[-which(nam==0)]}else{nums<-specnums}
  specprobs2[as.numeric(names(nums))]<-nums/sum(nums)
  #print(specprobs2)
  specprobs <- (1-params$poolprob)*specprobs2+params$poolprob*(1/params$nspec)
  #print(specprobs)
  specmat<-matrix(specprobs,nrow=params$nspec,ncol=params$ntype)
  pmat<-specmat*params$occmat
  #print(pmat)
  pmatvec<-pmat[,params$type[dead]]
  #print(pmatvec)
  y[dead]<-apply(pmatvec,2,function(x) which(rmultinom(1,1,x)==1)*rbinom(1,1,sum(x)))
  return(list(y))
}

params<-list()
params$nspec<-15
params$nsite<-1000
params$poolprob<-0.95
params$d<-0.1
params$occmat<-tmat
params$ntype<-12
params$type<-as.integer(gl(12,k=1000/12,length=1000))

test<-model2(1,rep(0,1000),params)



test<-ode(rep(0,1000),c(1:1000),model2,params,"iteration")

