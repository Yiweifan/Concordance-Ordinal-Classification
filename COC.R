norm2<-function(x){sqrt(sum(x^2))}

############linear learning############################
#objective function
fr<-function(t,train_X,train_Y,Xmatrix,k,Lambda,w,mu=0){
  
  n=nrow(train_X); p=ncol(train_X)
  
  beta=t[1:p]-t[-c(1:p)]
  beta_x=beta%*%t(train_X)
  
  temp1=beta_x[,train_Y<=k]
  temp2=beta_x[,train_Y>=(k+1)]
  temp3=c(outer(temp2,temp1,`-`))*sqrt(n)/10
  
  z=-sum(pnorm(temp3))/(n*(n-1))+Lambda*sum(w*t)+mu*(norm(beta,type="2"))^2
  
  return(z)
}

#gradient function
gr=function(t,train_X,train_Y,Xmatrix,k,Lambda,w,mu=0){
  
  n=nrow(train_X); p=ncol(train_X)
  
  beta=t[1:p]-t[-c(1:p)]
  beta_x=beta%*%t(train_X)
  
  temp3=beta_x[,train_Y<=k]
  temp4=beta_x[,train_Y>=(k+1)]
  temp5=dnorm(outer(temp4,temp3,`-`)*sqrt(n)/10)
  temp5=lapply(seq_len(nrow(temp5)),function(x){temp5[x,]})

  temp0=mapply("%*%",temp5,Xmatrix,SIMPLIFY=F)
  temp0=matrix(unlist(temp0),ncol=p,byrow=T)
  temp=-colSums(temp0)
  temp_final=rbind(temp+2*mu*beta,
                   -temp-2*mu*beta)

  z=(c(t(temp_final))*sqrt(n)/10)/(n*(n-1))+Lambda*w
  
  return(as.vector(z))
}

#maximize concordance function
solu<-function(initial,train_X,train_Y,k,Lambda,w,no.it=10,mu=0){
  p=ncol(train_X)
  
  temp1=train_X[train_Y<=k,]
  temp2=train_X[train_Y>=(k+1),]
  temp2=lapply(seq_len(nrow(temp2)),function(x){temp2[x,]})
  temp2=lapply(temp2,function(x){t(x-t(temp1))})
  
  t=optim(par=initial,fn=fr,gr=gr,
          train_X=train_X,train_Y=train_Y,
          Xmatrix=temp2,k=k,
          Lambda=Lambda,w=w,method="L-BFGS-B",mu=mu,
          lower=0,control=list(maxit=no.it))$par
  
  beta=t[1:p]-t[-c(1:p)]
  
  return(beta)
}


#evaluate classification
evaluation<-function(predy,testy,cumprob,K){
  #input:
    #predy: predicted labels
    #testy: true labels
    #cumprob: predicted cumulative probabilities
    #K: number of categories
  #output: MER, MAE, Kendall's coefficient, MRPS
  predy=as.numeric(predy); testy=as.numeric(testy)
  n=length(predy)
  
  c=0; d=0; et=0; ep=0 
  for(i in 1:n){
    c=c+sum(((predy-predy[i])*(testy-testy[i]))>0)
    d=d+sum(((predy-predy[i])*(testy-testy[i]))<0)
    et=et+sum((predy-predy[i])==0&(testy-testy[i])!=0)
    ep=ep+sum((testy-testy[i])==0&(predy-predy[i])!=0)
  }
  
  loss0_1=sum(predy!=testy)/length(testy)
  loss_ab=sum(abs(predy-testy))/length(testy)
  taub=(c-d)/(sqrt(c+d+et)*sqrt(c+d+ep))
  
  rps=lapply(seq_len(K-1),function(k){
    sum((cumprob[k,]-(testy<=k))^2)
    })
  loss_rps=sum(unlist(rps))/length(testy)
  return(c(loss0_1,loss_ab,1-taub,loss_rps))
}


#predict labels
mycut<-function(beta,test_X,cutoff){
  #input:
    #beta: estimated coefficients
    #test_X: covariates
    #cutoff: estimated thresholds
  #output: predicted labels
  p=ncol(test_X)
  score=test_X%*%t(beta)
  temp=t(t(score)-cutoff)
  temp=cbind(temp,-1)
  predy=apply(temp,1,function(x){min(which(x<0))})
  
  temp=temp>0
  temp=apply(temp,1,function(x){all(diff(x)<=0)})
  prop=1-sum(temp)/length(temp)
  return(list(predy=predy,prop=prop))
}


#predict cumulative probabilities
pred_prob<-function(beta,gamma,trainx,trainy,testx,K){
  #input:
    #beta: estimated coefficients
    #gamma: estimated thresholds
  #output: predicted cumulative probabilities
  
  score=t(gamma-t(trainx%*%t(beta)))
  tempy=c()
  for(k in 1:(K-1)){
    tempy=cbind(tempy,trainy>k)
  }
  lower=c(score);  lower=lower[c(tempy)!=0]
  upper=c(score);  upper=upper[c(tempy)!=1]
  lower=c(lower[which(lower>=min(upper))],
          max(lower[which(lower<min(upper))]))
  upper=c(upper[which(upper<=max(lower))],
          min(upper[which(upper>max(lower))]))
  lower=c(lower,rep(Inf,length(upper)))
  upper=c(rep(Inf,length(lower)-length(upper)),upper)
  constant=min(c(lower,upper))
  lower=lower-constant
  upper=upper-constant
  lower[lower==Inf]=0
  
  iserror=tryCatch({fit=ic_np(cbind(lower,upper))
  Ghatfun<-stepfun((fit$T_bull_Intervals[1,]+constant),
                   c(0,cumsum(fit$p_hat)),right=T)
  score=t(gamma-t(testx%*%t(beta)))
  cumprob=rbind(apply(score,1,Ghatfun),1)},
  error=function(e) e)
  
  if(inherits(iserror,"error")){
    cumprob=matrix(1,nrow=K,ncol=nrow(testx))
  }
  
  return(cumprob)
}


#minimize 0-1 loss
optimspliti_01<-function(beta,train_X,train_Y,K){
  #input:
    #beta: estimated coefficients
  #output: optimal thresholds
  n=nrow(train_X); p=ncol(train_X)
  
  myspliti=rep(0,(K-1))
  score=c(train_X%*%beta[1,])
  tempy=train_Y[order(score)]; tempy[tempy>1]=2
  score=sort(score)
  temp=matrix(1,n,n)
  temp[lower.tri(temp,diag=F)]=2
  err=apply(temp,2,function(x){sum(x!=tempy)/n})
  myspliti[1]=score[which.min(err)]
  predy=rep(1,which.min(err))
  
  if(length(predy)!=n){
    for(k in 2:(K-1)){
      score=c(train_X%*%beta[k,])
      tempy=train_Y[order(score)]; tempy[tempy>k]=k+1
      score=sort(score)
      temp=matrix(k,n,n)
      temp[lower.tri(temp,diag=F)]=(k+1)
      temp=temp[,-c(1:length(predy))]
      if(length(predy)==(n-1)){
        temp[1:length(predy)]=predy
        err=sum(temp!=tempy)/n
      }else{
        temp[1:length(predy),]=predy
        err=apply(temp,2,function(x){sum(x!=tempy)/n})
      }
      myspliti[k]=score[(which.min(err)+length(predy))]
      predy=c(predy,rep(k,sum(as.matrix(temp)[,which.min(err)]==k) ))
      if(length(predy)==n){break}
    }
  }
  
  return(myspliti)
}


#minimize abs loss
optimspliti_abs<-function(beta,train_X,train_Y,K){
  #input:
    #beta: estimated coefficients
  #output: optimal thresholds
  n=nrow(train_X); p=ncol(train_X)
  
  myspliti=rep(0,(K-1))
  score=c(train_X%*%beta[1,])
  tempy=train_Y[order(score)]; tempy[tempy>1]=2
  score=sort(score)
  temp=matrix(1,n,n)
  temp[lower.tri(temp,diag=F)]=2
  err=apply(temp,2,function(x){sum(abs(x-tempy))/n})
  myspliti[1]=score[which.min(err)]
  predy=rep(1,which.min(err))
  
  if(length(predy)!=n){
    for(k in 2:(K-1)){
      score=c(train_X%*%beta[k,])
      tempy=train_Y[order(score)]; tempy[tempy>k]=k+1
      score=sort(score)
      temp=matrix(k,n,n)
      temp[lower.tri(temp,diag=F)]=(k+1)
      temp=temp[,-c(1:length(predy))]
      if(length(predy)==(n-1)){
        temp[1:length(predy)]=predy
        err=sum(abs(temp-tempy))/n
      }else{
        temp[1:length(predy),]=predy
        err=apply(temp,2,function(x){sum(abs(x-tempy))/n})
      }
      myspliti[k]=score[(which.min(err)+length(predy))]
      predy=c(predy,rep(k,sum(as.matrix(temp)[,which.min(err)]==k) ))
      if(length(predy)==n){break}
    }
  }
  
  return(myspliti)
}


#minimize rps
optimspliti_rps<-function(initial,beta,train_X,train_Y,K){
  #input:
    #beta: estimated coefficients
  #output: optimal thresholds
  n=nrow(train_X); p=ncol(train_X)
  
  score=train_X%*%t(beta)
  tempscore=score
  myspliti=initial
  score=apply(score,2,sort)
  constrain=apply(as.matrix(score[,c(2:(K-1))]-score[,c(1:(K-2))]),2,max)
  tempspliti=initial
  tempspliti=matrix(rep(tempspliti,n),ncol=n)
  tempspliti[1,]=score[,1]
  err=apply(tempspliti,2,function(gamma){
    cumprob=pred_prob(beta,gamma,train_X,train_Y,train_X,K)
    rps=lapply(seq_len(K-1),function(k){
      sum((cumprob[k,]-(train_Y<=k))^2)
    })
    sum(unlist(rps))/length(train_Y)
  })
  myspliti[1]=score[which.min(err),1]
  tempspliti[1,]=myspliti[1]
  
  for(k in 2:(K-1)){
    check=((tempscore[,(k-1)]-myspliti[(k-1)])<=0)
    if(sum(check)==0){
      thre=Inf
    }else{thre=max(tempscore[which(check),k])}
    index=(score[,k]>=thre)
    #index=((score[,k]-myspliti[(k-1)])>=constrain[(k-1)])
    if(sum(index)==0){break}
    tempspliti[k,]=score[,k]
    if(sum(index)==1){
      err=0
    }else{
      err=apply(tempspliti[,index],2,function(gamma){
        cumprob=pred_prob(beta,gamma,train_X,train_Y,train_X,K)
        rps=lapply(seq_len(K-1),function(k){
          sum((cumprob[k,]-(train_Y<=k))^2)
        })
        sum(unlist(rps))/length(train_Y)
      })
    }
    myspliti[k]=score[(n-sum(index)+which.min(err)),k]
    tempspliti[k,]=myspliti[k]
  }
  
  return(myspliti)
}


#cross validation of COC
cv.cocnet<-function(trainx,trainy,K,weights,
                    lambda,nfolds,foldid,ini,no.it=10){
  p=ncol(trainx);  Cnum=length(lambda)
  
  index_list=cbind(rep(1:nfolds,each=Cnum),rep(1:Cnum,nfolds))
  cl=makeCluster(8)
  registerDoParallel(cl)
  cv_error=foreach(ii=1:nrow(index_list),.combine=c,.export=c("fr","gr","solu","norm2","optimspliti_01","mycut")) %dopar% {
    s=index_list[ii,1];    j=index_list[ii,2]
    
    rank_CI=lapply(seq_len(K-1),function(k){
      solu(c(ini[k,],rep(0,p)),trainx[-foldid[[s]],],trainy[-foldid[[s]]],
           k,lambda[j],weights[(1+2*(k-1)*p):(2*k*p)],no.it=no.it)
    })
    rank_CI=matrix(unlist(rank_CI),ncol=p,byrow=T)
    
    l2=apply(rank_CI,1,norm2);  l2[which(l2==0)]=1
    beta=rank_CI/l2
    c=optimspliti_01(beta,trainx[-foldid[[s]],],trainy[-foldid[[s]]],K)
    predy=mycut(beta,trainx[foldid[[s]],],c)$predy
    sum(predy!=trainy[foldid[[s]]])/length(predy)
  }
  stopCluster(cl)
  
  cv_error=matrix(cv_error,nrow=nfolds,byrow=T)
  
  return(list(lambda=lambda,cvm=cv_error))
}



#COC given the value of the tuning parameter
coc<-function(trainx,trainy,K,weights,lambda,
              loss="all",rescale=1,ini,no.it=10){
  p=ncol(trainx)
  
  cl=makeCluster(min((K-1),8))
  registerDoParallel(cl)
  rank_CI=foreach(k=1:(K-1),.combine=rbind,.export=c("fr","gr","solu")) %dopar% {
    solu(c(ini[k,],rep(0,p)),trainx,trainy,k,lambda,weights[(1+2*(k-1)*p):(2*k*p)],no.it=no.it)
  }
  stopCluster(cl)
  
  l2=apply(rank_CI,1,norm2);  l2[which(l2==0)]=1
  beta=rank_CI/l2*rescale
  
  gamma=list()
  if(loss=="01"){
    gamma[["01"]]=optimspliti_01(beta,trainx,trainy,K)
  }
  if(loss=="abs"){
    gamma[["abs"]]=optimspliti_abs(beta,trainx,trainy,K)
  }
  if(loss=="rps"){
    initial=optimspliti_abs(beta,trainx,trainy,K)
    gamma[["rps"]]=optimspliti_rps(initial,beta,trainx,trainy,K)
  }
  if(loss=="all"){
    gamma[["01"]]=optimspliti_01(beta,trainx,trainy,K)
    gamma[["abs"]]=optimspliti_abs(beta,trainx,trainy,K)
    gamma[["rps"]]=optimspliti_rps(gamma[["abs"]],beta,trainx,trainy,K)
  }
  
  return(list(coefficients=beta,thresholds=gamma))
}


