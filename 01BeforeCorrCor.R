####association before correction
rawdata=readRDS("C:\\lxx\\metha\\Meth2022\\MethRaw.rds")
rawdata=rawdata[sort(rownames(rawdata)),]
phen=read.table("C:\\lxx\\metha\\Meth2022\\Age-Sex.txt",header=T,sep="\t",stringsAsFactors = F)
Mfemal=phen[phen$SexDetailed%in%c("female","female-young"),"GSM"]
Mmal=phen[phen$SexDetailed%in%c("Male","male"),"GSM"]
Mcas=phen[phen$SexDetailed%in%c("male-castrated"),"GSM"]
Mmale=apply(rawdata[rownames(rawdata),c(Mmal,Mcas)],1,mean)
Mfemale=apply(rawdata[rownames(rawdata),Mfemal],1,mean)
#Mcastrated=apply(rawdata[rownames(rawdata)%in%uniq,Mcas],1,mean)

setwd("C:\\lxx\\metha\\Meth2022\\PositiveControls")
f=c(70,75,80,85,90,95)
corresult=as.data.frame(f)
samp=unlist(strsplit(grep("F70cgid.txt",list.files(),value = T),".F70cgid.txt"))
for(s in 1:length(samp)){
for(i in 1:length(f)){
  name=paste0(samp[s],".F",f[i],"cgid.txt")
  data=read.table(name,header=T,sep="\t",stringsAsFactors = F)
  rownames(data)=data$CGid
  data=data[sort(data$CGid),]
  #data1=data[rownames(data)%in%clockid,]
  cp=c(1:10,15,20,25,30,35,40,45,50)
  fval=paste0("S",cp,"_F")
  for(n in 1:length(fval)){
    if(samp[s]%in%ancmale){
      corresult[i,fval[n]]=cor(Mmale[rownames(data)],data[,fval[n]],method="spearman")
    }else{
      corresult[i,fval[n]]=cor(Mfemale[rownames(data)],data[,fval[n]],method="spearman")
    }
  }
  colnames(corresult)=c("f",fval)
}

write.table(corresult,paste0(samp[s],"corr.before.txt"),col.names = T,row.names = F,sep="\t",quote=F)
pdf(paste0(samp[s],"corr.before.pdf"))
plotdata=read.table(paste0("GVA602_AMIS_1_00818corr.before.txt"),header=T,sep="\t",stringsAsFactors = F)
cp=c(1:10,15,20,25,30,35,40,45,50)
fval=paste0("S",cp,"_F")
col=rcolors::get_color(rcolors::rcolors$t2m_29lev,n=length(cp))
col=colorRampPalette(c("#AFEEEE","navyblue"))(length(cp))
pdf("GVA602_AMIS_1_00818-SF50-wgs.pdf")
plot(plotdata$f,plotdata$S1_F,type="o",pch=18,ylim=c(0.4,0.9),col=col[1],xlab="Celler methylation fraction percent",ylab="R")
for(m in 2:length(cp)){
  lines(plotdata$f,plotdata[,fval[m]],type="o",pch=18,col=col[m])
}
lab=gsub("F","CpG",gsub("S","",colnames(plotdata)[-1]))
legend("top",fval,inset=-0.01,lab,col=col,pch=rep(19,length(cp)),lwd=2,bty="n",ncol=5,cex=1)
dev.off()
}