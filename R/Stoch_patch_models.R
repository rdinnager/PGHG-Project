stoch_patch<-function(t,y,params) {
  empt<-y==0 #determine empty patches, 0 means empty, any other y means species with ID number 'y'
  #Add new mortality to empty patch vector; d=death rate
  dead<-c(sample(which(!empt),rbinom(1,sum(!empt),params$d)),which(empt)) 
  specnums<-table(y[!empt]) #count how many there are of each species
  specprobs2<-rep(0,params$nspec) #initialize species probability vector
  nam<-as.numeric(names(specnums)) #grab the ID numbers of species
  #if (0 %in% nam){nums<-specnums[-which(nam==0)]}else{nums<-specnums} #Get rid of empty cells
  specprobs2[as.numeric(names(nums))]<-nums/sum(nums) #Species prob is equal to proportion of cells occupied
  #add equal probability of being drawn from species pool. poolprob = prob of immigrant coming from pool
  specprobs <- (1-params$poolprob)*specprobs2+params$poolprob*(1/params$nspec) 
  specmat<-matrix(specprobs,nrow=params$nspec,ncol=params$ntype) #but probs into species*resource matrix
  pmat<-specmat*params$occmat #multiple prob of immigrating by prob of being able to feed (0 or 1)
  pmatvec<-pmat[,params$type[dead]] #pull out probs for all cells according to their resource type
  #Use probs to calculate which species should be placed into which empty cells if any
  y[dead]<-apply(pmatvec,2,function(x) which(rmultinom(1,1,x)==1)*rbinom(1,1,sum(x)))
  return(list(y))
}