library(sest)
rawdata=readRDS("C:\\lxx\\metha\\Meth2022\\MethRaw.rds")
rawdata=rawdata[order(rownames(rawdata)),]
chrX=read.table("C:\\lxx\\metha\\Meth2022\\PositiveControls\\chrX.txt",header=T,sep="\t",stringsAsFactors = F);chrX=chrX[,1]
datall=rawdata[chrX,];dim(datall)
data=t(datall)
phen=read.table("C:\\lxx\\metha\\Meth2022\\Age-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
mmale=phen[phen$SexDetailed%in%c("male","male-castrated"),"GSM"]
mfemale=phen[phen$SexDetailed%in%c("female"),"GSM"]
data219=(readRDS("C:\\lxx\\metha\\Meth2022\\PositiveControls\\68Merged-dis0-dep50-lmall.RDS"))
data=(data219[1:68,colnames(data219)%in%chrX])
sex102=read.table("C:/lxx/metha/Meth2022/PositiveControls/88Indv-Sex.txt",header=T,sep="\t",stringsAsFactors = F)

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
datall=data219[names(bb1[bb1<0.1])[names(bb1[bb1<0.1])%in%chrX],names(bb[bb<0.1])];dim(datall)
datalll=cbind(datall,rawdata[rownames(datall),])

library(impute)
imputelob1= impute.knn(as.matrix(datalll),k = 2, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069) 
imputedlob1 <- as.matrix(imputelob1$data)
data=t(imputedlob1)
data=data[rownames(data)%in%paste0("predic",sex102$oldsample),]
sex102=sex102[sex102$oldsample%in%gsub("predic","",rownames(data)) & sex102$DataType=='WGS',]
sex102=sex102[order(sex102$oldsample),]

phendata=phen[,c("SexDetailed","GSM")]
phendata=sex102[,c("Sex","oldsample")]
colnames(phendata)=c("sex","sample")
phendata$sample=gsub("predic","",phendata$sample)
Xestsex=function(data,intervals,phendata,return_with_reference=FALSE){
rownames(data)=gsub("predic","",rownames(data))
pca.all.X <- prcomp((merge_intervals(data,intervals)), scale.=FALSE)
kmeans.X <- kmeans(pca.all.X$x, 2)
reference.male=phendata[phendata$sex%in%c("male","male-castrated","M"),"sample"]#paste0("predic",phen[phen$Sex%in%"M","sample"])#
reference.female=phendata[phendata$sex%in%c("female","F"),"sample"]#paste0("predic",phen[phen$Sex%in%"F","sample"])#
reference.cas=phendata[phendata$sex%in%"male-castrated","GSM"]
samples.all=rownames(data)
cluster_name.male.X <- names(table(kmeans.X$cluster[reference.male]))[order(table(kmeans.X$cluster[reference.male]), decreasing=TRUE)[1]]
cluster_name.female.X <- names(table(kmeans.X$cluster[reference.female]))[order(table(kmeans.X$cluster[reference.female]), decreasing=TRUE)[1]]
tab.pred <- data.frame(row.names=samples.all, X.PC1=pca.all.X$x[samples.all, "PC1"], X.PC2=pca.all.X$x[samples.all, "PC2"],predicted.X=NA)
tab.pred[, "predicted.X"] <- ifelse(kmeans.X$cluster[samples.all] == cluster_name.male.X, "M", "F")
tab.pred[, "true"] <- unlist(lapply(1:length(samples.all),function(x){phendata[phendata$sample==samples.all[x],"sex"]}))
sex_estimation <- list(test=tab.pred[samples.all, c("X.PC1","predicted.X")] )

if (return_with_reference) {
  sex_estimation$reference <- tab.pred[samples.all, c("X.PC1","X.PC2","predicted.X","true")]
}
tab.pred
}

merge_intervals <- function(data = NULL, intervals=NULL){
 tab.prop=lapply(1:nrow(data),function(m){
  tab.base_prop=as.data.frame((data[m,]))
  tab.prop <- data.frame(row.names=rownames(data)[m])
  prev <- as.numeric(intervals[1])
  for (i in 2:length(intervals)) {
    curr <- as.numeric(intervals[i])
    interval_name <- sprintf("(%g,%g]", prev, curr)
    tab.prop[, interval_name] <- length(tab.base_prop[tab.base_prop>prev &tab.base_prop<=curr])/ncol(data)
    prev <- curr
  }
  return(tab.prop)
  })
  do.call(rbind,tab.prop)
}
mergeintvnew=function(data,intervals){
  rownames(data)=gsub("predic","",rownames(data))
  tab.prob=do.call(rbind, lapply(rownames(data), function(x) table(cut(data[x,], breaks=intervals, include.lowest=FALSE))/ncol(data)))
  rownames(tab.prob)=rownames(data)
  tab.prob
}
library(reshape2)
library(ggplot2)
plotDist=function(data,intervals,phendata,filename){
  reference.male=phendata[phendata$sex%in%c("M","male","male-castrated"),"sample"]#phen[phen$SexDetailed%in%c("male","male-castrated"),"GSM"]
  reference.female=phendata[phendata$sex%in%c("F","female"),"sample"]
datamerge=merge_intervals(data,intervals)
rownames(datamerge)=gsub("predic","",rownames(datamerge))
g.ref.X.beta <- ggplot() + ylim(0,0.55) + theme_bw()
datamerge.melted <- melt(cbind(sample=rownames(datamerge), datamerge))
colnames(datamerge.melted)=c("sample","range","proportion")
datamerge.melted$sample=as.factor(datamerge.melted$sample)
datamerge.melted$range=factor(datamerge.melted$range,levels=colnames(datamerge))
datamerge.melted$sex=factor(ifelse(datamerge.melted$sample %in% reference.male, "Male", "Female"), levels=c("Male", "Female"))
g.ref.X.beta <- g.ref.X.beta + geom_line(data=datamerge.melted, aes(x=range, y=proportion, group=sample, col=sex), alpha=0.02, linetype="solid") +
  scale_color_discrete(type=c("#000066","#F8766D"))+
  geom_boxplot(data=datamerge.melted, aes(x=range, y=proportion, fill=sex), alpha=0.2)+
   scale_fill_manual(values=c("Male" = "#000066", "Female" = "#F8766D"))+
  theme(legend.position = c(0.2,0.9),legend.key.size=unit(0.3,"inches"))
ggsave(paste0(filename,"-chrX-sexdistribution.pdf"),width=200,height=200,units="mm")
}
datawrong=tab.pred[c("GSM5326985","GSM5326987"),]
plotPCA=function(tab.pred,filename){
  library(ggplot2)
  ggplot(tab.pred,aes(as.numeric(X.PC1),as.numeric(X.PC2),color=true,alpha=0.8))+geom_point(size=3)+
    xlab(paste0("PC1(",round(sd(tab.pred$X.PC1),4)*100,"%)"))+ylab(paste0("PC2(",round(sd(tab.pred$X.PC2),4)*100,"%)"))+
    scale_color_manual(values=c("#F8766D","#8B4500","#000066","#FFCC00","#556B2F","blue"))
  ggsave(filename,width=200,height=197,units="mm")
}

##pca for all
dat=t(imputedlob1)
data=dat[c(paste0("predic",sex102$oldsample),phen$GSM),]
rownames(data)=c(paste0("predic",sex102$oldsample),phen$GSM)
phendata=phen[,c("SexDetailed","GSM")]
colnames(phendata)=c("sex","sample")
phendata1=sex102[,c("Sex","oldsample")]
colnames(phendata1)=c("sex","sample")
phendata=rbind(phendata1,phendata)
tab=Xestimatesex(data,seq(0,1,0.2),phendata,return_with_reference = T)
tab$true=gsub("F","ancF",tab$true)
tab$true=gsub("M","ancM",tab$true)
wrong=tab[tab$true%in%c("male","male-castrated") & tab$X.PC1<0.05,]
tab=rbind(tab,wrong,wrong,wrong)
ggplot(tab,aes(as.numeric(X.PC1),as.numeric(X.PC2),color=true,alpha=0.8))+geom_point(size=3)+
  xlab(paste0("PC1(",round(sd(tab$X.PC1),4)*100,"%)"))+ylab(paste0("PC2(",round(sd(tab$X.PC2),4)*100,"%)"))+
  scale_color_manual(values=c("#8B4500","#556B2F","#F8766D","#000066","#FFCC00","blue"))+theme_bw()
ggsave("chrX-mdanc-960.pdf",width=200,height=197,units="mm")
