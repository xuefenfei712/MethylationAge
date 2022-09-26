data219=(readRDS("C:\\lxx\\metha\\Meth2022\\PositiveControls\\68Merged-dis0-dep50-lmall.RDS"))
chrX=read.table("chrX.txt",header=T,sep="\t",stringsAsFactors = F);chrX=chrX[,1];length(chrX)
rawdata=readRDS("C:\\lxx\\metha\\Meth2022\\MethRaw.rds")
rawdata=rawdata[order(rownames(rawdata)),]
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
datall=cbind(data219[names(bb1[bb1<0.1]),names(bb[bb<0.1])],rawdata[names(bb1[bb1<0.1]),gsm]);dim(datall)
library(impute)
imputelob1= impute.knn(as.matrix(datall),k = 2, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069) 
imputedlob1 <- as.matrix(imputelob1$data)
datall=imputedlob1
#datall=data219[rownames(data219)%in%chrX==FALSE,names(bb[bb<0.1])];dim(datall)

phen=read.table("C:\\lxx\\metha\\Meth2022\\Age-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
gsm=phen$GSM#[phen$Tissue%in%c("Blood","Liver"),"GSM"]
datallX=datall[rownames(datall)%in%chrX,];dim(datallX)
datallnoX=datall[rownames(datall)%in%chrX==FALSE,];dim(datallnoX)
sex102=read.table("C:/lxx/metha/Meth2022/PositiveControls/88Indv-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
sex102=sex102[sex102$DataType%in%"WGS",];dim(sex102)

datallnoXx=datallnoX[,colnames(datallnoX)%in%c(paste0("predic",sex102$oldsample),gsm)];dim(datallnoXx)
datallXx=datallX[,colnames(datallX)%in%c(paste0("predic",sex102$oldsample),gsm)];dim(datallXx)
tab=data.frame("")
tab$sex=NULL
tab$ratio=NULL
for(i in 1:length(grep("predic",colnames(datallnoXx)))){
  tab[i,"sex"]=sex102[sex102$oldsample%in%gsub("predic","",colnames(datallnoXx)[i]),"Sex"]
  tab[i,"ratio"]=mean(datallnoXx[,i],na.rm=T)/mean(datallXx[,i],na.rm=T)
  }
for(i in (length(grep("predic",colnames(datallnoXx)))+1):ncol(datallnoXx)){
  tab[i,"sex"]=phen[phen$GSM%in%(colnames(datallnoXx)[i]),"SexDetailed"]
  tab[i,"ratio"]=mean(datallnoXx[,i],na.rm=T)/mean(datallXx[,i],na.rm=T)
}
tab=tab[,2:3];head(tab)
rownames(tab)=colnames(datallnoXx)
tab$sex=gsub("M","ancM",tab$sex)
tab$sex=gsub("F","ancF",tab$sex)
tab$tis=c(rep("anc",50),phen$Tissue)
tab$cat=paste0(tab$sex,"-",tab$tis)
write.table(tab,"tab.txt",col.names = T,row.names = T,sep='\t',quote=F)
tab=read.table("tab.txt",header = T,sep="\t",stringsAsFactors = F)
library(ggsignif)
library(ggplot2)
p = ggplot(tab, aes(x=cat2, y=ratio,color=sex)) + 
  geom_boxplot()+#scale_x_discrete(limits=c("ancF","ancM","female","male","male-castrated"))+
  theme_bw()+theme(axis.text.x=element_text(angle=75, hjust=1),axis.text.y=element_text(size=5,hjust=0))+
  ylab("Autosome-to-X ratio")+geom_signif(comparisons = tab$sex,step_increase = 0.1,map_signif_level = F,test = t.test)
ggsave("sexRatio50-960-333-Bliver.pdf",width=200,height=200,units="mm")
