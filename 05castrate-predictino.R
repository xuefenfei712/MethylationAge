library("WGCNA")
library("GO.db")
library(reshape2)
library(stringr)
rawdata=readRDS("C:\\lxx\\metha\\Meth2022\\MethRaw.rds")
rawdata=rawdata[sort(rownames(rawdata)),]
age=read.table("C:\\lxx\\metha\\Meth2022\\Age-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
age[age$GSM%in%c("GSM5326985","GSM5326987"),"SexDetailed"]=rep("female",2)
m=c("male","male-castrated")
casage=11
Mm=age[age$SexDetailed%in%m &age$Tissue%in%c("Blood") & age$Age>casage ,"GSM"]
Fm=age[age$SexDetailed%in%"female" &age$Tissue%in%c("Blood") &age$Age>casage,"GSM"]
TdataM=t(rawdata[,colnames(rawdata)%in%Mm])
TdataF=t(rawdata[,colnames(rawdata)%in%Fm])
ageM=age[age$SexDetailed%in%m &age$Tissue%in%c("Blood") &age$Age>casage,]
ageF=age[age$SexDetailed%in%"female" &age$Tissue%in%c("Blood") &age$Age>casage,]

sex=c(ageM$SexDetailed)
options(stringsAsFactors = T)# ??????????????????,??????????????????,?????????????????????,????????????????????????FALSE
enableWGCNAThreads()#???????????????
datExpr <- TdataM
weight=ifelse(sex%in%"male-castrated",2,1)
sexresult=standardScreeningBinaryTrait(datExpr,weight,kruskalTest = T,qValues = T,na.action = "na.exclude",getAreaUnderROC = T)
sigMC=sexresult[abs(sexresult$corPearson)>=0.5,"ID"];length(sigMC)
######ancient data
merg13=readRDS("C:\\lxx\\metha\\Meth2022\\PositiveControls\\68Merged-dis0-dep50-lmall-64-29113.RDS")

sex102=read.table("C:/lxx/metha/Meth2022/PositiveControls/102Indv-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
age34=read.table("68Merged-dis0-dep50-lmall.143Bliver-CI-wgs-linear.txt",header=T,sep="\t",stringsAsFactors = F)
age34=age34[age34$Tissue%in%"anc",]
ancm=paste0("predic",age34[(age34$sex%in%c("M","F") & age34$DNAmAgeloo>11),"name"])
ancmm=paste0("predic",age34[(age34$sex%in%c("M") & age34$DNAmAgeloo>11),"name"])
dataA=t(merg13[,ancm])
sexA=age34[paste0("predic",age34$name)%in%rownames(dataA),"sex"]
datal=(cbind(rawdata[rownames(merg13),],merg13))
sexall=c(ageM$SexDetailed,ageF$SexDetailed,sexA)

library(impute)
imputelob1= impute.knn(as.matrix((datal)),k = 2, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069) 
imputedlob1 <- as.matrix(imputelob1$data)
datall=t(imputedlob1)
datallm=datall[c(Mm,Fm,ancm),]
pca.all<- prcomp(datall[c(Mm,Fm,ancm),colnames(datall)], scale.=FALSE)
pcatab <- kmeans(pca.all$x, 4)
pcatab=as.data.frame(cbind(pca.all$x[,1:2],sexall,c(ageM$Age,ageF$Age,age34[paste0("predic",age34$name)%in%ancm,"DNAmAgeloo"])))
colnames(pcatab)=c("PC1","PC2","sexall","age")
library(ggplot2)
pcall=ggplot(pcatab,aes(as.numeric(PC1),as.numeric(PC2),color=as.numeric(age),shape=sexall))+geom_point()+
  xlab(paste0("PC1(",round(sd(pcatab$PC1),2),"%)"))+ylab(paste0("PC2(",round(sd(pcatab$PC2),2),"%)"))+
  geom_polygon(data = plyr::ddply(pcatab, "sexall", function(df) df[chull(df[[1]], df[[2]]), ]), fill = NA)+
  scale_color_gradient2(low = "#F4A582",mid="blue",high = "darkred")+ theme_bw()
ggsave("PCA-sex-Blood-MA11.pdf",width=200,height=200,units="mm")

set.seed(100)
cas=datall[c(Mm),colnames(datall)%in%sigMC]
cassex=ageM$SexDetailed
nmds <- vegan::metaMDS(cas, distance = 'bray',trace=2)
result <- nmds$points
result <- as.data.frame(cbind(result, rownames(result)))
#?????? https://www.sci666.com.cn/39320.html
colnames(result)[1:3] <- c("NMDS1", "NMDS2", "samples")
result$NMDS1 <- as.numeric(as.character(result$NMDS1))
result$NMDS2 <- as.numeric(as.character(result$NMDS2))
result$samples <- as.character(result$samples)
result <- cbind(result, c(ageM$SexDetailed),c(ageM$Age))
colnames(result) <- c("NMDS1", "NMDS2", "samples","sex","Age")
head(result)
library(ggplot2)
#berpred=result[27:34,]
mc=ggplot(result, aes(NMDS1, NMDS2, colour=Age,shape=sex)) +
  geom_polygon(data = plyr::ddply(result, "sex", function(df) df[chull(df[[1]], df[[2]]), ]), fill = NA) +
  geom_point(size=5)+scale_color_gradient(low = "blue",high = "darkred")+
  scale_shape_manual(values=c("-","+"))+ theme_bw()
 # annotate("text",x=berpred$NMDS1+0.01,y=berpred$NMDS2-0.01,label=berpred$samples,col="red",size=3)
ggsave(paste0(length(sigMC),"-MDS-RF-sex-Blood-NEW11-0.55.pdf"),width=200,height=200,units="mm")

####male less than 11
Mm11=age[age$SexDetailed%in%m &age$Tissue%in%c("Blood") &age$Age,"GSM"]
Mm11sex=age[age$SexDetailed%in%m &age$Tissue%in%c("Blood")&age$Age,"SexDetailed"]
Mm11age=age[age$SexDetailed%in%m &age$Tissue%in%c("Blood")&age$Age,"Age"]

cas11=datall[Mm11,colnames(datall)%in%sigMC]
nmds <- vegan::metaMDS(cas11, distance = 'bray',na.rm=T)
result <- nmds$points
result <- as.data.frame(cbind(result, rownames(result)))
#?????? https://www.sci666.com.cn/39320.html
colnames(result)[1:3] <- c("NMDS1", "NMDS2", "samples")
result$NMDS1 <- as.numeric(as.character(result$NMDS1))
result$NMDS2 <- as.numeric(as.character(result$NMDS2))
result$samples <- as.character(result$samples)
result <- cbind(result, Mm11sex,Mm11age)
colnames(result) <- c("NMDS1", "NMDS2", "samples","sex","Age")
head(result)
library(ggplot2)
#berpred=result[27:34,]
mc11=ggplot(result, aes(NMDS1, NMDS2, colour=Age,shape=sex)) +
  geom_polygon(data = plyr::ddply(result, "sex", function(df) df[chull(df[[1]], df[[2]]), ]), fill = NA) +
  geom_point(size=5)+scale_color_gradient2(low = "#F4A582",mid = "blue",high="darkred")+
  scale_shape_manual(values=c("-","+"))+ theme_bw()
# annotate("text",x=berpred$NMDS1+0.01,y=berpred$NMDS2-0.01,label=berpred$samples,col="red",size=3)
ggsave(paste0(length(sigMC),"-MDS-RF-sex-Blood-les11-0.55.pdf"),width=200,height=200,units="mm")


set.seed(123)
RFsex <- randomForest::randomForest(datall[Mm,colnames(datall)%in%sigMC],as.factor(ageM$SexDetailed), cv.fold=10,ntree=200, importance = TRUE)
RFsex

#get best mtry and ntree
n=ncol(cas)-1
rate=NULL
for(i in 1:n){
  set.seed(123)
  RFsex1 <- randomForest::randomForest(datall[Mm,colnames(datall)%in%sigMC],as.factor(ageM$SexDetailed), mtry=i,ntree=200, proximity=TRUE,importance = TRUE)
  rate[i]=mean(RFsex1$err.rate,na.rm = T)
}
mtry=which.min(rate)
set.seed(123)
RFsex2<- randomForest::randomForest(datall[Mm,colnames(datall)%in%sigMC],as.factor(ageM$SexDetailed), mtry=mtry,ntree=60, proximity=TRUE,importance = TRUE)
plot(RFsex2,col=1:1)
predic=predict(RFsex2,datall[,colnames(datall)%in%sigMC])
plot(randomForest::margin(RFsex2, predic), main = "truth rate")

a=as.data.frame(randomForest::importance(RFsex2))
sigmcsig=rownames(a[a$MeanDecreaseAccuracy>=1,])

ind=sample(2,length(Mm),replace=TRUE,prob=c(0.8,0.2))
datall1=datall[c(Mm),colnames(datall)%in%sigMC]
datallc=datall[c(ancmm),colnames(datall)%in%sigMC]
train=datall1[ind==1,]
test=datall1[ind==2,]
trainsex=ageM[ind==1,"SexDetailed"]
testsex=ageM[ind==2,"SexDetailed"]
set.seed(1)
RFsex3<- randomForest::randomForest(datall1,as.factor(ageM$SexDetailed), mtry=mtry,ntree=20, proximity=TRUE,importance = TRUE)
rfsex3predic=predict(RFsex3,test)
library("pROC")
pred<-prediction(rfsex3predic[,2],testsex)
perf<-performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)
multiclass.roc(as.ordered(testsex),as.ordered(rfsex3predic))
roc<-roc(as.ordered(testsex) ,as.ordered(rfsex3predic))
plot(roc,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)

pred1rf <- predict(RFsex3,datallc)
pred2rf <- predict(RFsex3,rbind(datall1,datallc), type="prob")

tab=as.data.frame(pred2rf)
tab1=cbind(gsub("predic","",rownames(tab)),tab,rep("rf",nrow(tab)))
colnames(tab1)=c("name","male","male-castrated")
write.table(tab1,"rf-predictSex0.6.txt",col.names = T,row.names = F,sep="\t",quote=F)

#SVM
library(e1071)
mod <- svm(datall1,as.factor(ageM$SexDetailed), kernel = c("linear"),
           cost = 2^(-9),type="C",gamma = 1,cross=10,probability = TRUE,ntree = 500)

suppVec=rownames(mod$SV)
pred1 <- predict(mod,test)
pred2 <- predict(mod,test, probability = TRUE)
pred1 <- predict(mod,datallc)
pred2 <- predict(mod,rbind(datall1,datallc), probability = TRUE)
res <- as.matrix(attr(pred2, "probabilities")[, c("male-castrated","male")])
ress=cbind(rownames(res),res,rep("SVM",nrow(res)))
write.table(ress,"SVM-predic-Sex0.6.txt",col.names = T,row.names = F,sep="\t",quote=F)

###gknn
modknn=gknn(datall1,as.factor(ageM$SexDetailed),k=sqrt(nrow(train)),method="manhattan")
knnclas=predict(modknn, test, type = c("class"))
knnraw=predict(modknn,test, type = c("prob"))
knnclas=predict(modknn, rbind(datall1,datallc), type = c("class"))
knnraw=predict(modknn,rbind(datall1,datallc), type = c("prob"))
knntab=cbind(rownames(knnraw),knnraw,rep("gknn",nrow(knnraw)))
colnames(knntab)=c("name","male","male-castrated","method")
write.table(knntab,"knn-predict-Sex0.6.txt",col.names = T,row.names = F,sep="\t",quote=F)

###annotate
anndata=read.csv("C:\\lxx\\metha\\DesignFiles\\Equus_caballus.equcab3.0.100.HorvathMammalMethylChip40.v1.csv",header=T,sep=",",stringsAsFactors = F)
siggene=anndata[anndata$CGid%in%colnames(cas),c("CGid","SYMBOL","GeneRegionID", "seqnames","CGstart" ,"CGend")]
write.table(siggene,"sig122-castrated0.5.txt",col.names = T,row.names = T,sep="\t")
