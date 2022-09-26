library(doMC)
registerDoMC()

rftest=function(x){
  names(datamerge)=c("SF","SDep","SDis","Mfemale")
    testid=sort(sample(nrow(datamerge),round(nrow(datamerge)*0.1),0))
    testdata=datamerge[testid,]
    traindata=datamerge[-testid,]
    set.seed(123)
    rfmodel=foreach(ntree=rep(50,10), .combine=randomForest::combine,
              .multicombine=TRUE, .packages='randomForest') %dopar% {
 randomForest::randomForest(I(Mfemale^x) ~ (SF)^2+SF+SDep+SF*SDep,type= "regression",cv.fold=10,data=traindata,ntree=ntree,importance=TRUE)
}
    aft=cbind(datamerge$Mfemale,predict(rfmodel,datamerge))
    colnames(aft)=c("Mfemale","predict")
 aftd=as.data.frame(aft[abs(aft[,1]-aft[,2])<sd(aft[,1]),])#residuals(lm(Mfemale~predict,as.data.frame(aft)))),])
   MX <- density(aftd[,1])$x
   MY <- density(aftd[,1])$y
   Mpeak <- MX[which(diff(sign(diff(MY)))==-2)]
   AX <- density(aftd[,2])$x
    AY <- density(aftd[,2])$y
    Apeak <- MX[which(diff(sign(diff(AY)))==-2)]
  abs(min(Apeak)-min(Mpeak))+abs(max(Apeak)-max(Mpeak))-cor(aftd[,1],aftd[,2],method="spearman")#+abs(max(Apeak)-min(Apeak)-max(Mpeak)+min(Mpeak))
    }

###R value
rfmd=function(datamerge,x){
  if(nrow(datamerge>0)){
    testid=sample(1:nrow(datamerge),round(nrow(datamerge)*0.2),0)
    testdata=datamerge[testid,]
    traindata=datamerge[-testid,]
    rfmodel=randomForest::randomForest((Mfemale)^x~(SF)^2+SDep,cv.fold=10,data=traindata,ntree=500,importance=TRUE)
    cor(datamerge$Mfemale,predict(rfmodel,datamerge))
  }else{
    0
  }
}
##predict matrix
rfpred=function(datamerge,dep,x){
  datamerge=datamerge[datamerge$SDep>dep,]
  names(datamerge)=c("SF","SDep","SDis","Mfemale")
  if(nrow(datamerge>0)){
    testid=sample(1:nrow(datamerge),round(nrow(datamerge)*0.1),0)
    testdata=datamerge[testid,]
    traindata=datamerge[-testid,]
 rfmodel=foreach(ntree=rep(50,10), .combine=randomForest::combine,
              .multicombine=TRUE, .packages='randomForest') %dopar% {
 randomForest::randomForest(I(Mfemale^x) ~ (SF)^2+SF+SDep+SF*SDep,type= "regression",cv.fold=10,data=datamerge,ntree=ntree,importance=TRUE)
}
    aft=cbind(datamerge$Mfemale,predict(rfmodel,datamerge))
    colnames(aft)=c("Mfemale","predict")
resmd=lm(Mfemale~predict,as.data.frame(aft))
    aa=cbind(aft[,1],predict(resmd))
    aftd=aa[abs(aa[,1]-aa[,2])<sd(aa[,1]),]
    colnames(aftd)=c("Modern","predict")
    as.data.frame(aftd[,2])
  }else{
    0
  }
}
rfmdplot=function(datamerge,dep,x){
  names(datamerge)=c("SF","SDep","SDis","Mfemale")
  datamerge=datamerge[datamerge$SDep>dep,]
  if(nrow(datamerge>0)){
    trainid=sample(seq(1:nrow(datamerge)),round(nrow(datamerge)*0.9,0))
    traindata=datamerge[trainid,]
    testdata=datamerge[-trainid,]
 rfmodel=foreach(ntree=rep(50,10), .combine=randomForest::combine,
              .multicombine=TRUE, .packages='randomForest') %dopar% {
 randomForest::randomForest(I(Mfemale^x) ~ (SF)^2+SF+SDep+SF*SDep,type= "regression",cv.fold=10,data=datamerge,ntree=ntree,importance=TRUE)
}
    aft=cbind(datamerge$Mfemale,predict(rfmodel,datamerge))
    colnames(aft)=c("Modern","predict")
 resmd=lm(Modern~predict,as.data.frame(aft))
    aa=cbind(aft[,1],predict(resmd))
    aftd=aa[abs(aa[,1]-aa[,2])<sd(aft[,1]),]  
      colnames(aftd)=c("Modern","predict")           
    psych::pairs.panels(as.data.frame(aftd),method="spearman",breaks=30,main=paste("AftCor",dep,nrow(aftd),sep="-"))
  }
}

