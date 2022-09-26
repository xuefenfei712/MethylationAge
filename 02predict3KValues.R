clock=read.table("C:\\lxx\\metha\\Meth2022\\clockCGid.txt",header=F,sep="\t",stringsAsFactors = F)
clockid=sort(clock[,1]);length(clockid)
sigid=read.table("C:/lxx/metha/Meth2022/PositiveControls/2022nc-age-sex-sigCGid.txt",header=T,sep="\t",stringsAsFactors = F)
uniq=read.table("C:\\lxx\\metha\\Meth2022\\uniqclock1.txt",header=F,sep="\t")
uniq=uniq[,1];length(uniq)
rawdata=readRDS("C:\\lxx\\metha\\Meth2022\\MethRaw.rds")
rawdata=rawdata[sort(rownames(rawdata)),]
phen=read.table("C:\\lxx\\metha\\Meth2022\\Age-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
sex102=read.table("C:/lxx/metha/Meth2022/PositiveControls/102Indv-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
ancm=sex102[sex102$Sex%in%"M","sample"]
ancf=sex102[sex102$Sex%in%"F","sample"]
set.seed(100)
Mfemal=phen[phen$SexDetailed%in%c("female","female-young"),"GSM"]
Mfemale=apply(rawdata[rownames(rawdata),Mfemal],1,mean)
Mmal=phen[phen$SexDetailed%in%c("Male","male"),"GSM"]
Mcas=phen[phen$SexDetailed%in%c("male-castrated"),"GSM"]
Mmale=apply(rawdata[rownames(rawdata),c(Mmal,Mcas)],1,mean)
args <- commandArgs(T)
best=read.table("C:/lxx/metha/Meth2022/PositiveControls/bestSF.txt",header=F,sep="\t",stringsAsFactors = F)
source("rfcorrection-para.R")
for(b in 1:nrow(best)){
  data90=read.table(paste0(best[b,1],".F95cgid.txt"),header=T,sep="\t",stringsAsFactors = F)
  rownames(data90)=data90$CGid
  sf=best[b,2]
  sd=gsub("F","FDepth",best[b,2])
  sdis=gsub("F","DistFCpG",best[b,2])
  data1=data90[rownames(data90),c(sf,sd,sdis)]
   data1=data1[sort(rownames(data1)),]
 if(best[b,1]%in%ancm){
   Mmale1=Mmale[rownames(data1)]
  datald=cbind(data1,Mmale1)
  colnames(datald)=c("SF","SDep","SDis","Mmale")
  }else{
    Mfemale1=Mfemale[rownames(data1)]
    datald=cbind(data1,Mfemale1)
    colnames(datald)=c("SF","SDep","SDis","Mfemale")
  }
  ld=corldplot(datald,assoc=best[b,1])
  data1=data1[as.numeric(data1[,sdis])<=0,]
 # data1=data1[sort(rownames(data1)),]
 if(best[b,1]%in%ancm){
   Mmale1=Mmale[rownames(data1)]
  datamerge=cbind(data1,Mmale1);
  datamerge=datamerge[rownames(datamerge)%in%sigid[,1]==FALSE,];print(dim(datamerge))
  colnames(datamerge)=c("SF","SDep","SDis","Mmale")
  }else{
    Mfemale1=Mfemale[rownames(data1)]
    datamerge=cbind(data1,Mfemale1);
    datamerge=datamerge[rownames(datamerge)%in%sigid[,1]==FALSE,];print(dim(datamerge))
    colnames(datamerge)=c("SF","SDep","SDis","Mfemale")
  }
 para=lapply(1:10,function(i){
  pa=optimize(rftest,c(-3,3),maximum = F)
  para=pa$minimum
  })
  x=mean(as.numeric(para),na.rm=T)
  pred=(rfpred(datamerge,dep=100,x=x))
if(best[b,1]%in%ancm){
  datamerge1=cbind(data1,Mmale[rownames(data1)]);dim(datamerge1)
  colnames(datamerge1)=c("SF","SDep","SDis","Mmale")
  }else{
    datamerge1=cbind(data1,Mfemale[rownames(data1)]);dim(datamerge1)
    colnames(datamerge1)=c("SF","SDep","SDis","Mfemale")
  }
 #pred=(rfpred(datamerge1,dep=50,x=x))
  names(pred)=paste0("predic",best[b,1])
saveRDS(pred,paste0(best[b,1],"pred31k-dis0-dep50.RDS"))

if(best[b,1]%in%ancm){
data=datamerge1[,c("Mmale","SF")]
}else{
  data=datamerge1[,c("Mfemale","SF")] 
}
colnames(data)=c("Modern","SF")
pdf(paste0(best[b,1],"-AftCorDis0Dep50",nrow(data),".pdf"))
psych::pairs.panels(as.data.frame(data),method="spearman",breaks=30,main=paste0(best[b,1],"_M95_",best[b,2],"-",nrow(data)))
#psych::pairs.panels(data,main=paste(best[b,1],"AftCor",nrow(aa),sep="-"),method="spe",breaks=30)
rfmdplot(datamerge1,50,x)
rfmdplot(datamerge1,100,x)
rfmdplot(datamerge1,150,x)
rfmdplot(datamerge1,200,x)
#psych::pairs.panels(aa,main=paste(best[b,1],"AftCor",nrow(aa),sep="-"))
print(paste0(best[b,1]," is done!"))
dev.off()
}
