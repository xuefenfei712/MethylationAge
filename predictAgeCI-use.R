args <- commandArgs(T)
data219=readRDS(args[1])
data219=data219[,grep("predic",colnames(data219))]
#data219=data2191[-which(apply(data2191[,1:length(grep("predic",colnames(data2191)))],1,function(x) all(is.na(x)))),]
rawdata=readRDS("/disk/regine/data2/xuexue/Methalization/predicValue/MethRaw.rds")
rawdata=rawdata[sort(rownames(rawdata)),]

calNA=function(data){
a=NULL;b=NULL
for(i in 1:ncol(data)){
 a=length(which(is.na(data[,i])))
          b=c(b,a)
}
names(b)=colnames(data)
  b
}
bb=calNA(data219)/nrow(data219)
calsiteNA=function(data){
  a=NULL;b=NULL
  for(i in 1:nrow(data)){ 
    a=length(which(is.na(data[i,])))
    b=c(b,a)
  }
  names(b)=rownames(data)
  b
}
bb1=calsiteNA(data219)/ncol(data219)

datall=cbind(data219[(names(bb1[bb1<0.1])),names(bb[bb<0.1])],rawdata[names(bb1[bb1<0.1]),]);dim(datall)
library(impute)##col is individual, row is site
imputelob1= impute.knn(as.matrix((datall)),k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069) 
rawdclock<- as.matrix(t(imputelob1$data))
rawdclockk=t(apply(rawdclock,1,function(x) scale(x,center = T,scale = T)))
colnames(rawdclockk)=colnames(rawdclock)
###scale
rownames(rawdclockk)=rownames(rawdclock)
#colnames(rawdclockk)=colnames(rawdclock)
rawdclock=(rawdclockk);dim(rawdclock)

print(paste0("0. there are sites:",ncol(rawdclock)))

saveRDS(rawdclock[grep("predic",rownames(rawdclock),value = T),],gsub("RDS","imp.RDS",args[1]))
###############clock
tis=c("Blood","Liver")
phen=read.table("/disk/regine/data2/xuexue/Methalization/predicValue/Age-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
age=phen[phen$Tissue%in%tis,]
data=rawdclockk[age$GSM,]

agespan=read.table("GvawithAgespan.txt",header=T,sep="\t",stringsAsFactors = F)
rownames(agespan)=agespan$sample
agespan1=agespan[agespan$DataType%in%"WGS" & agespan$clockin%in%"yes",];dim(agespan)
agespan
ageall=c(agespan1$age,age$Age);length(agespan1$age)
all=rbind(rawdclock[paste0("predic",agespan1$sample),],data);dim(all)
#print age
a=as.data.frame(cbind(grep("predic",rownames(rawdclock),value = T),rep(0,length(grep("predic",rownames(rawdclock))))))
for(i in 1:nrow(a)){
  a[i,2]=ifelse(any(agespan$sample%in%gsub("predic","",a[i,1])),agespan[agespan$sample%in%gsub("predic","",a[i,1]),"age"],0)
  }
colnames(a)=c("ID","age")
write.table(a,paste0(nrow(rawdclock),"ageraw.txt"),col.names=T,row.names=F,sep=" ",quote=F)

library(glmnet)
library(doMC)
registerDoMC()
set.seed(10)
cvfit = glmnet::cv.glmnet(all, ageall, family="gaussian", type.measure = "mse",standardize = FALSE,alpha=0.9,keep=T,nfolds=10,parallel=TRUE,trace.it = TRUE)
bestlambda=cvfit$lambda.min
fit=glmnet(all,ageall, family="gaussian", type.measure = "mse",standardize = FALSE,alpha=0.9,keep=T,nfolds=10,trace.it = TRUE,nlambda=100)
cor=as.matrix(predict(fit, type="coef", s=bestlambda))
sig=names(cor[cor[,1]!=0,]);length(sig)
cortab=as.data.frame(cbind(names(cor[cor[,1]!=0,]),cor[cor[,1]!=0,]));dim(cortab)
colnames(cortab)=c("probe","coef")
#write.table(cortab,paste0(nrow(rawdclock),"indv-bliver-cap.coef"),col.names=T,row.names=F,sep=" ",quote=F)

rawdclock=rawdclock[,colnames(rawdclock)%in%sig]
age=phen[phen$Tissue%in%tis,]
data=rawdclock[age$GSM,]
ageuse=agespan[agespan$DataType%in%"WGS","sample"]
all1=rbind(rawdclock[paste0("predic",ageuse),],(data))
age1=c(agespan[ageuse,"age"],age$Age)

library(foreach)
lmsample=function(all,i){
set.seed(i)
 (sort(sample(seq(ncol(all)),ncol(all),replace = T)))
}

ela=NULL
elas=matrix("")
elasm=matrix("",ncol=3)
for(i in 1:nrow(all1)){
  ela=lapply(1:100,function(m){
    set.seed(i+m)
    ss=lmsample(all1,i+m)
    cvfit = glmnet::cv.glmnet(all1[-i,ss], age1[-i], family="gaussian", type.measure = "mse",standardize = FALSE,alpha=0.9,keep=T,nfolds=10,parallel=TRUE,trace.it = TRUE)
    bestlambda=cvfit$lambda.min
    fit1=glmnet(all1[-i,ss],age1[-i], family="gaussian", type.measure = "mse",standardize = FALSE,alpha=0.9,keep=T,nfolds=10,trace.it = TRUE,nlambda=100)
    ela=(predict(fit1, all1[i,ss], type="link", s=bestlambda))
  })
  elas=t(data.frame(Rmisc::CI(unlist(ela),ci=0.95)))
  rownames(elas)=rownames(all1)[i]
  elasm=rbind(elasm,elas)
}
elasm=elasm[-1,]
sex=c(agespan[ageuse,"sex"],age$SexDetailed)
tis=c(rep("anc",length(ageuse)),age$Tissue)
 type=c(rep("anc",length(ageuse)),rep("md",nrow(age)))
 result=as.data.frame(cbind(rownames(elasm),age1,sex,elasm,tis,type))
 colnames(result)=c("sample","Age","Sex","upper","DNAmAgeloo","lower","Tissue","type")
 result$Age=as.numeric(result$Age)
 result$DNAmAgeloo=as.numeric(result$DNAmAgeloo)
 result1=result[result$type%in%"anc",]
  for(i in 1:nrow(result1)){
  result1[i,"agest"]=agespan[agespan$sample%in%gsub("predic","",result1[i,"sample"]),"agest"]
  result1[i,"ageend"]=agespan[agespan$sample%in%gsub("predic","",result1[i,"sample"]),"ageend"]
  }
mae=round(mean(abs(as.numeric(result[result$type%in%"anc","Age"])-as.numeric(result[result$type%in%"anc","DNAmAgeloo"]))),2)
cor=round(cor(as.numeric(result[result$type%in%"anc","Age"]),as.numeric(result[result$type%in%"anc","DNAmAgeloo"])),2)
print(paste0(mae,cor))
library(ggplot2)
p=ggplot(result,aes(x=Age,y=DNAmAgeloo,colour=Tissue))+geom_point()+
  geom_abline(intercept = 0,slope=1,linetype="dashed")+#scale_color_manual(values=c("#FF0000","#CD853F","#556B2F"))+
  geom_smooth(method="lm",se=F,size=0.6,color="black",aes(group=1))+
  annotate("segment",x=result1$agest,xend=result1$ageend,y=(result1$DNAmAge),yend=(result1$DNAmAge))+
  annotate("segment",x=as.numeric(result1$Age),xend=as.numeric(result1$Age),y=as.numeric(result1$lower),yend=as.numeric(result1$upper),col="grey")+
  geom_text(aes(x=5,y=25,label=paste0("cor=",cor," mae=",mae)),parse=F,show.legend=NA)+
  scale_fill_identity(guide = "legend")
ggsave(paste0(gsub("RDS","",args[1]),length(sig)-1,"BLiver-CI-cap-linear.pdf"),width=200,height=200,units="mm")
write.table(result,paste0(gsub("RDS","",args[1]),length(sig)-1,"Bliver-CI-cap-linear.txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(result1,paste0(gsub("RDS","",args[1]),length(sig)-1,"Bliver-CI-cap-linear-agespan.txt"),col.names=T,row.names=F,sep="\t",quote=F)
#################predict unknown
unk=rawdclock[rownames(rawdclock)%in%agespan$sample==FALSE,]
unkdata=unk[rownames(unk)%in%grep("predic",rownames(unk),value = T),colnames(unk)%in%sig];dim(unkdata)
print(paste0("there are unknown sites: ",nrow(unkdata)))
all3=all1[,colnames(all1)%in%sig];dim(all3)
###clock
elu=NULL
elus=matrix("")
elusm=matrix("",ncol=3)
for(i in 1:nrow(unkdata)){
  elu=lapply(1:length(sig),function(m){
    set.seed(i+m)
    ss=lmsample(all3,i+m)
cvfit = glmnet::cv.glmnet(all3[,ss], age1, family="gaussian", type.measure = "mse",standardize = FALSE,alpha=0.9,keep=T,nfolds=10,parallel=TRUE,trace.it = TRUE)
bestlambda=cvfit$lambda.min
fit=glmnet(all3[,ss],age1, family="gaussian", type.measure = "mse",standardize = FALSE,alpha=0.9,keep=T,nfolds=10,trace.it = TRUE,nlambda=100)
elu=(predict(fit, unkdata[i,ss], type="link", s=bestlambda))
})
  elus=t(data.frame(Rmisc::CI(unlist(elu),ci=0.95)))
  rownames(elus)=rownames(unkdata)[i]
  elusm=rbind(elusm,elus)
}
elusm=elusm[-1,]
sex102=read.table("/disk/regine/data2/xuexue/Methalization/88Indv-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
sex1021=sex102[order(sex102$sample),]
d=cbind(grep("predic",rownames(unk),value = T),sex1021[sex1021$sample%in%gsub("predic","",rownames(unkdata)),"Sex"],elusm,rep("anc",nrow(unkdata)))
write.table(as.data.frame(d),paste0(gsub("RDS","",args[1]),length(sig),"BLiverunknown-cap.txt"),col.names=T,row.names=F,sep="\t",quote=F)

