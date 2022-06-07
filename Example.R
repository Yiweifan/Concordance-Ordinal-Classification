#load packages
library(caret)
library(foreach)
library(doParallel)
library(icenReg)
library(dummies)
source("COC.R")
set.seed(345)


#load data
mydata=read.table("balance.txt",sep=",")
y=mydata$V1;  y=match(y,c("L","B","R"))
x=mydata[,-1];  x=scale(x)

K=max(y); p=ncol(x); n=nrow(x)
mydata=cbind(x,y); mydata=as.matrix(mydata)
colnames(mydata)=c(paste("x",c(1:p),sep=""),"y")

#split training and test sets
train_row=sample(c(1:n),(n/3))
trainx=mydata[train_row,1:p]
trainy=mydata[train_row,(p+1)]

testx=mydata[-c(train_row),1:p]
testy=mydata[-c(train_row),(p+1)] 

#creat 5-folds
folds=createFolds(c(1:nrow(trainx)),k=5,list=TRUE,returnTrain=FALSE)

#regularization path 
Cnum=41; Cvec=10^seq(-7,-3,length.out=Cnum)
#adaptive lasso weights
#setup parallel backend to use many processors
cl=makeCluster(detectCores()-1) #not to overload your computer
registerDoParallel(cl)
w=foreach(k=1:(K-1),.combine=rbind,.export=c("fr","gr","solu")) %dopar% {
  solu(rep(0,(2*p)),trainx,trainy,k,Lambda=0,rep(1,(2*p)),no.it=100,mu=10^(-6))
}
stopCluster(cl)
w=t(apply(w,1,function(x){x/norm2(x)}));  w=1/abs(w)
w=w[rep(seq_len(nrow(w)),each=2),];  w=c(t(w))
w[w==Inf]=10^4

#5-folds cross validation
cocfit=cv.cocnet(trainx,trainy,K,w,lambda=Cvec,
                 nfolds=5,foldid=folds,ini=matrix(0,nrow=(K-1),ncol=p),no.it=100)
cv_error=cocfit$cvm
#best lambda: one standard error rule
temp=colMeans(cv_error);  mysd=apply(cv_error,2,sd)
onese=min(temp)+mysd[which.min(temp)]/sqrt(5)
best_lambda=Cvec[max(which(temp<onese))]

#initialize cumulative probabilities
cumprob=matrix(1,nrow=K,ncol=length(testy))
#train on the whole training set with three different losses
cocfit=coc(trainx,trainy,K,w,best_lambda,loss="all",
           ini=matrix(0,nrow=(K-1),ncol=p),no.it=100)
beta=cocfit$coefficients

#0-1 loss
gamma=cocfit$thresholds[["01"]]
predy=mycut(beta,testx,gamma)$predy
errcoc=evaluation(predy,testy,cumprob,K)[1]

#abs loss
gamma=cocfit$thresholds[["abs"]]
predy=mycut(beta,testx,gamma)$predy
errcoc=c(errcoc,evaluation(predy,testy,cumprob,K)[2:3])

#rps
gamma=cocfit$thresholds[["rps"]]
cumprob=pred_prob(beta,gamma,trainx,trainy,testx,K)
errcoc=c(errcoc,evaluation(predy,testy,cumprob,K)[4])

#evaluation metric: MER, MAE, Kendall's coefficient, MRPS
print(errcoc)
