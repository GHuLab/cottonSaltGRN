########################################叶（时间+品种）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample =fread("C:/Users/hr345/Desktop/all/E各类sample/sample(final).txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6]) #batch:批次
dim(expr)
#[1] 75376   170
anno=sample
dim(anno)
#234 6

#####叶   去掉明显离群样本
anno=anno[-140,]
anno=anno[-197:-202,]
anno=anno[-1:-123,]
dim(anno)
#104 6

expr1=t(expr)
expr1=expr1[-140,]
expr1=expr1[-197:-202,]
expr1=expr1[-1:-123,]
expr=t(expr1)
dim(expr)

keep=rowSums(expr>1) >= floor(0.05*ncol(expr))
table(keep)
#keep
#FALSE  TRUE 
#20586 54790
expr <- expr[keep,]
dim(expr)
expr1=t(expr)

#PCA
library(ggplot2)
df_pca <- prcomp(log2(expr1+1)) #计算主成分
dim(expr1)
# 162 75376
anno$time<- factor(anno$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession)  
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  # xlab(percentage[1]) +
  # ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))

#tSNE
library(Rtsne)
tsne_out = Rtsne(expr1,perplexity=20)
pdat = data.frame(tsne_out$Y)
pdat=cbind(pdat,anno$time,anno$Accession)
colnames(pdat) = c("tSNE_1","tSNE_2","time","Accession")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
pdat$Accession<- factor(pdat$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))

#UMAP
library(umap)
iris.umap = umap::umap(expr1)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$time,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","time","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h","1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))
									   
###去批次
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
dim(ad.batch)

ad.factors.df <- data.frame(batch = anno$batch)
expr1=expr1
class(expr1) 
ad.clr=expr1
#install.packages("vegan")
library(vegan)
library(limma)
library(mixOmics)
#removeBatchEffect
ad.mod <- model.matrix( ~ anno$batch)
dim(ad.clr)
ad.rBE <- t(removeBatchEffect(t(ad.clr), batch = anno$batch,design = ad.mod))

#PCA
library(gridExtra)
library("ggpubr")
library(ggplot2)
df_pca <- prcomp(ad.rBE) #计算主成分
dim(ad.rBE)
#104 75376
dim(expr1)
#104 75376
anno$time<- factor(anno$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession)  

head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))

#tSNE
library(Rtsne)
tsne_out = Rtsne(ad.rBE,perplexity=20)
pdat = data.frame(tsne_out$Y)
pdat=cbind(pdat,anno$time,anno$Accession)
colnames(pdat) = c("tSNE_1","tSNE_2","time","Accession")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
pdat$Accession<- factor(pdat$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))

#UMAP
library(umap)
iris.umap = umap::umap(ad.rBE)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$time,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","time","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h","1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))


#####第二种方法
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("sva")
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))

#PCA
df_pca <- prcomp(ad.ComBat) #计算主成分
dim(ad.ComBat)
#104 75376
dim(expr1)
#104 75376
anno$time<- factor(anno$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession)  

head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))

#tSNE
library(Rtsne)
tsne_out = Rtsne(ad.ComBat,perplexity=20)
pdat = data.frame(tsne_out$Y)
pdat=cbind(pdat,anno$time,anno$Accession)
colnames(pdat) = c("tSNE_1","tSNE_2","time","Accession")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
pdat$Accession<- factor(pdat$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))
  
#UMAP
library(umap)
iris.umap = umap::umap(ad.ComBat)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$time,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","time","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h","1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))


###########叶（批次）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/E各类sample/sample(final).txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6]) #batch:批次
dim(expr)
#[1] 75376   170
anno=sample
dim(anno)
#234 6

#####叶   去掉明显离群样本
anno=anno[-140,]
anno=anno[-197:-202,]
anno=anno[-1:-123,]
dim(anno)
#104 6

expr1=t(expr)
expr1=expr1[-140,]
expr1=expr1[-197:-202,]
expr1=expr1[-1:-123,]
expr=t(expr1)
dim(expr)

keep=rowSums(expr>1) >= floor(0.05*ncol(expr))
table(keep)
#keep
#FALSE  TRUE 
#20483 54893
expr <- expr[keep,]
dim(expr)
expr1=t(expr)
##未去批次
library(ggplot2)
df_pca <- prcomp(log2(expr1+1)) #计算主成分
dim(expr1)
# 162 75376
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)  
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+
  geom_point(size=3)+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#fc913a","#ff4e50","#cbe86b","#3b8686"))


#去批次
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
dim(ad.batch)
expr1=expr1
class(expr1) 
ad.clr=expr1
library(vegan)
library(limma)
library(mixOmics)

#removeBatchEffect
ad.mod <- model.matrix( ~ anno$batch)
dim(ad.clr)
ad.rBE <- t(removeBatchEffect(t(ad.clr), batch = anno$batch,design = ad.mod))
library(gridExtra)
library("ggpubr")
library(ggplot2)
df_pca <- prcomp(ad.rBE) #计算主成分
dim(ad.rBE)
#104 75376
dim(expr1)
#104 75376
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)  
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+ geom_point()
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+
  geom_point(size=3)+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#fc913a","#ff4e50","#cbe86b","#3b8686"))


#combat
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))
df_pca <- prcomp(ad.ComBat) #计算主成分
dim(ad.ComBat)
#104 75376
dim(expr1)
#104 75376
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)  

head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+
  geom_point(size=3)+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#fc913a","#ff4e50","#cbe86b","#3b8686"))






###########根（时间+品种）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
library(xlsx)
sample <- read.xlsx("C:/Users/hr345/Desktop/all/E各类sample/sample（去批次）.xlsx",1,encoding="UTF-8")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6]) #batch:批次
dim(expr)
#[1] 75376   170

anno=sample
dim(anno)
anno=anno[1:123,]
dim(anno)

expr1=t(expr)
expr1=expr1[1:123,]
dim(expr1)
expr1=t(expr1)
dim(expr1)
#[1]  75376 123
keep=rowSums(expr1>1) >= floor(0.05*ncol(expr1))
table(keep)
#FALSE  TRUE 
#18613 56763
expr <- expr1[keep,]
dim(expr)
expr1=t(expr)

#PCA
library(ggplot2)
df_pca <- prcomp(log2(expr1+1)) #计算主成分
dim(expr1)
anno$time<- factor(anno$time, levels=c("0h", "15min", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession) 
head(df_pcs,3)
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+ geom_point()
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+
     geom_point(size=3)+ 
     xlab(percentage[1]) +
     ylab(percentage[2])+
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
          legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#88dba3","#a5dff9","#1ec0ff","#0080ff","#4ea1d3","#2b90d9"))

#tSNE
library(Rtsne)
tsne_out = Rtsne(expr1,perplexity=20)
pdat = data.frame(tsne_out$Y)
pdat=cbind(pdat,anno$time,anno$Accession)
colnames(pdat) = c("tSNE_1","tSNE_2","time","Accession")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
pdat$Accession<- factor(pdat$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#77919d","#5e5e5f","#88dba3","#60c5ba","#cbe86b","#4f953b","#1ec0ff","#0080ff","#4ea1d3","#eb9f9f","#bf209f","#6a60a9","#5c196b"))

#UMAP
library(umap)
iris.umap = umap::umap(expr1)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$time,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","time","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#77919d","#5e5e5f","#88dba3","#60c5ba","#cbe86b","#4f953b","#1ec0ff","#0080ff","#4ea1d3","#eb9f9f","#bf209f","#6a60a9","#5c196b"))

###去批次	 
anno$time<- factor(anno$time, levels=c("0h", "15min", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$batch<- factor(anno$batch, levels=c("1","8","9","11"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("1","8","9","11"), ordered=TRUE)
dim(ad.batch)
ad.trt = factor(anno$time, levels=c("0h", "15min", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
ad.factors.df <- data.frame(time = anno$time, batch = anno$batch)
expr1=expr1
class(expr1)
ad.clr=expr1

library(vegan)
library(limma)
library(mixOmics)
#removeBatchEffect
ad.mod <- model.matrix( ~ anno$batch)
dim(ad.clr)
ad.rBE <- t(removeBatchEffect(t(ad.clr), batch = anno$batch,design = ad.mod))

#PCA
library(gridExtra)
library("ggpubr")
library(ggplot2)
df_pca <- prcomp(ad.rBE) #计算主成分
dim(ad.rBE)
dim(expr1)
anno$time<- factor(anno$time, levels=c("0h", "15min", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession)  

head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+ geom_point()
head(df_pcs,3)
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+ geom_point()
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+
     geom_point(size=3)+ 
     xlab(percentage[1]) +
     ylab(percentage[2])+
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
          legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#88dba3","#a5dff9","#1ec0ff","#0080ff","#4ea1d3","#2b90d9"))

#tSNE
library(Rtsne)
tsne_out = Rtsne(ad.rBE,perplexity=20)
str(tsne_out)
pdat = data.frame(tsne_out$Y)
rownames(pdat)=anno$SRR
pdat=cbind(pdat,anno$time,anno$Accession)

colnames(pdat) = c("tSNE_1","tSNE_2","time","Accession")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
pdat$Accession<- factor(pdat$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#77919d","#5e5e5f","#88dba3","#60c5ba","#cbe86b","#4f953b","#1ec0ff","#0080ff","#4ea1d3","#eb9f9f","#bf209f","#6a60a9","#5c196b"))

#UMAP
library(umap)
iris.umap = umap::umap(ad.rBE)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)

pumap=cbind(pumap,anno$time,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","time","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#77919d","#5e5e5f","#88dba3","#60c5ba","#cbe86b","#4f953b","#1ec0ff","#0080ff","#4ea1d3","#eb9f9f","#bf209f","#6a60a9","#5c196b"))

#combat去批次
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))

#PCA
df_pca <- prcomp(ad.ComBat) #计算主成分
dim(ad.ComBat)
dim(expr1)
anno$time<- factor(anno$time, levels=c("0h", "15min", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession)
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=Accession))+
     geom_point(size=3)+ 
     xlab(percentage[1]) +
     ylab(percentage[2])+
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
           legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#88dba3","#a5dff9","#1ec0ff","#0080ff","#4ea1d3","#2b90d9"))

#tSNE	 
library(Rtsne)
tsne_out = Rtsne(ad.ComBat,perplexity=20)
str(tsne_out)
pdat = data.frame(tsne_out$Y)
rownames(pdat)=anno$SRR
pdat=cbind(pdat,anno$time,anno$Accession)

colnames(pdat) = c("tSNE_1","tSNE_2","time","Accession")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
pdat$Accession<- factor(pdat$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#77919d","#5e5e5f","#88dba3","#60c5ba","#cbe86b","#4f953b","#1ec0ff","#0080ff","#4ea1d3","#eb9f9f","#bf209f","#6a60a9","#5c196b"))

#UMAP
library(umap)
iris.umap = umap::umap(ad.ComBat)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)

pumap=cbind(pumap,anno$time,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","time","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h", "15min", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#88dba3","#a5dff9","#1ec0ff","#0080ff","#4ea1d3","#2b90d9"))

########################################根（浓度+品种）

rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
library(xlsx)
sample <- read.xlsx("C:/Users/hr345/Desktop/all/E各类sample/sample（去批次）.xlsx",1,encoding="UTF-8")
#sample <- read.xlsx("C:/Users/hr345/Desktop/all/E各类sample/sample（去批次）2.xlsx",1,encoding="UTF-8")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6]) #batch:批次
dim(expr)

anno=sample
dim(anno)
anno=anno[1:123,]
dim(anno)

expr1=t(expr)
expr1=expr1[1:123,]
dim(expr1)
expr1=t(expr1)
dim(expr1)
#[1]  75376 123
keep=rowSums(expr1>1) >= floor(0.05*ncol(expr1))
table(keep)
#FALSE  TRUE 
#18613 56763
expr <- expr1[keep,]
dim(expr)
expr1=t(expr)

anno$concentration<- factor(anno$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
anno$batch<- factor(anno$batch, levels=c("1","8","9","11"), ordered=TRUE)
#anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("1","8","9","11"), ordered=TRUE)
#ad.batch =  factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
dim(ad.batch)
ad.trt = factor(anno$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
ad.factors.df <- data.frame(concentration = anno$concentration, batch = anno$batch)
expr1=expr1
class(expr1) 
ad.clr=expr1

#install.packages("vegan")
library(vegan)
library(limma)
library(mixOmics)
#removeBatchEffect
ad.mod <- model.matrix( ~ anno$batch)
dim(ad.clr)
ad.rBE <- t(removeBatchEffect(t(ad.clr), batch = anno$batch,design = ad.mod))

#PCA
library(gridExtra)
library("ggpubr")
library(ggplot2)
df_pca <- prcomp(ad.rBE) #计算主成分
dim(ad.rBE)
dim(expr1)
anno$concentration<- factor(anno$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,concentration = anno$concentration,Accession=anno$Accession)
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=concentration,shape=Accession))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,color=concentration,shape=Accession))+
  geom_point(size=3)+ 
  # xlab(percentage[1]) +
  # ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#d62a9d","#6a60a9","#f199bc"))

#tSNE
library(Rtsne)
tsne_out = Rtsne(ad.rBE,perplexity=20)
str(tsne_out)
pdat = data.frame(tsne_out$Y)
rownames(pdat)=anno$SRR
pdat=cbind(pdat,anno$concentration,anno$Accession)

colnames(pdat) = c("tSNE_1","tSNE_2","concentration","Accession")
head(pdat)
pdat$concentration<- factor(pdat$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
pdat$Accession<- factor(pdat$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=concentration,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#d62a9d","#6a60a9","#f199bc"))

#UMAP
library(umap)
iris.umap = umap::umap(ad.rBE)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)

pumap=cbind(pumap,anno$concentration,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","concentration","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$concentration<- factor(pumap$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=concentration,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#d62a9d","#6a60a9","#f199bc"))

#combat去批次
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))

#PCA
df_pca <- prcomp(ad.ComBat) #计算主成分
dim(ad.ComBat)
dim(expr1)
anno$concentration<- factor(anno$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,concentration = anno$concentration,Accession=anno$Accession)
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=concentration,shape=Accession))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=concentration,shape=Accession))+
     geom_point(size=3)+ 
     # xlab(percentage[1]) +
     # ylab(percentage[2])+
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
           legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#d62a9d","#6a60a9","#f199bc"))

#tSNE
library(Rtsne)
tsne_out = Rtsne(ad.ComBat,perplexity=20)
str(tsne_out)
pdat = data.frame(tsne_out$Y)
rownames(pdat)=anno$SRR
pdat=cbind(pdat,anno$concentration,anno$Accession)

colnames(pdat) = c("tSNE_1","tSNE_2","concentration","Accession")
head(pdat)
pdat$concentration<- factor(pdat$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
pdat$Accession<- factor(pdat$Accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=concentration,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#d62a9d","#6a60a9","#f199bc"))
  
#UMAP
library(umap)
iris.umap = umap::umap(ad.ComBat)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)

pumap=cbind(pumap,anno$concentration,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","concentration","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$concentration<- factor(pumap$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
#pumap$batch<- factor(pumap$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=concentration,shape=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#d62a9d","#6a60a9","#f199bc")) 

#######根（时间+浓度）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/E各类sample/sample(final).txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6]) #batch:批次
dim(expr)
#[1] 75376   170
anno=sample
dim(anno)
#234 6
anno=sample
dim(anno)
anno=anno[1:123,]
dim(anno)

expr1=t(expr)
expr1=expr1[1:123,]
dim(expr1)
expr1=t(expr1)
dim(expr1)
#[1]  75376 123
keep=rowSums(expr1>1) >= floor(0.05*ncol(expr1))
table(keep)
#FALSE  TRUE 
#18613 56763
expr <- expr1[keep,]
dim(expr)
expr1=t(expr)

#PCA
library(ggplot2)
df_pca <- prcomp(log2(expr1+1)) #计算主成分
dim(expr1)
anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$concentration<- factor(anno$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,concentration=anno$concentration)
head(df_pcs,3)
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=concentration))+ geom_point()
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=concentration))+
    geom_point(size=3)+ 
    # xlab(percentage[1]) +
    # ylab(percentage[2])+
    theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
          legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
    theme(axis.title.y = element_text(size = 15))+ 
    scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#4ea1d3","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(15,14,17,18,16))
	
#tSNE
library(Rtsne)
tsne_out = Rtsne(expr1,perplexity=20)
pdat = data.frame(tsne_out$Y)
pdat=cbind(pdat,anno$time,anno$concentration)
colnames(pdat) = c("tSNE_1","tSNE_2","time","concentration")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
pdat$concentration<- factor(pdat$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=concentration))+
     geom_point(size=3)+ 
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
          legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#ef5285","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(15,14,17,18,16))

#UMAP
library(umap)
iris.umap = umap::umap(expr1)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$time,anno$concentration)
colnames(pumap) = c("UMAP_1","UMAP_2","time","concentration")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
pumap$concentration<- factor(pumap$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)
ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=concentration))+
     geom_point(size=3)+ 
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
           legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#ef5285","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(15,14,17,18,16))


#去批次
anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$concentration<- factor(anno$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
dim(ad.batch)
#ad.trt = factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
ad.factors.df <- data.frame(time = anno$time, batch = anno$batch, concentration = anno$concentration)
expr1=expr1
class(expr1)
ad.clr=expr1

##removeBatchEffect
library(vegan)
library(limma)
library(mixOmics)
#removeBatchEffect
ad.mod <- model.matrix( ~ anno$batch)
dim(ad.clr)
ad.rBE <- t(removeBatchEffect(t(ad.clr), batch = anno$batch,design = ad.mod))

#PCA
library(gridExtra)
library("ggpubr")
library(ggplot2)
df_pca <- prcomp(ad.rBE) #计算主成分
dim(ad.rBE)
dim(expr1)
anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$concentration<- factor(anno$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,concentration=anno$concentration)  

head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=concentration))+ geom_point()
head(df_pcs,3)
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=concentration))+ geom_point()
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=concentration))+
     geom_point(size=3)+ 
     xlab(percentage[1]) +
     ylab(percentage[2])+
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
          legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#4ea1d3","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(15,14,17,18,16))

#tSNE
library(Rtsne)
tsne_out = Rtsne(ad.rBE,perplexity=20)
str(tsne_out)
pdat = data.frame(tsne_out$Y)
rownames(pdat)=anno$SRR
pdat=cbind(pdat,anno$time,anno$concentration)

colnames(pdat) = c("tSNE_1","tSNE_2","time","concentration")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
pdat$concentration<- factor(pdat$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=concentration))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#ef5285","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(15,14,17,18,16))

#UMAP
library(umap)
iris.umap = umap::umap(ad.rBE)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)

pumap=cbind(pumap,anno$time,anno$concentration)
colnames(pumap) = c("UMAP_1","UMAP_2","time","concentration")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
pumap$concentration<- factor(pumap$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=concentration))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#ef5285","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(15,14,17,18,16))

#combat去批次
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))

#PCA
df_pca <- prcomp(ad.ComBat) #计算主成分
dim(ad.ComBat)
dim(expr1)
anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$concentration<- factor(anno$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,concentration=anno$concentration)
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=concentration))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=concentration))+
     geom_point(size=3)+ 
     xlab(percentage[1]) +
     ylab(percentage[2])+
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
           legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#4ea1d3","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(15,14,17,18,16))

#tSNE	 
library(Rtsne)
tsne_out = Rtsne(ad.ComBat,perplexity=20)
str(tsne_out)
pdat = data.frame(tsne_out$Y)
rownames(pdat)=anno$SRR
pdat=cbind(pdat,anno$time,anno$concentration)

colnames(pdat) = c("tSNE_1","tSNE_2","time","concentration")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
pdat$concentration<- factor(pdat$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=concentration))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#ef5285","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(15,14,17,18,16))

#UMAP
library(umap)
iris.umap = umap::umap(ad.ComBat)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)

pumap=cbind(pumap,anno$time,anno$concentration)
colnames(pumap) = c("UMAP_1","UMAP_2","time","concentration")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
pumap$concentration<- factor(pumap$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=concentration))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#D1B6E1","#fab1ce","#EE7785","#8EC0E4","#C5E99B","#8CD790","#77AF9C","#379392","#285943"))+ scale_shape_manual(values=c(18,17,15,14,16))

#######根（批次）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/E各类sample/sample(final).txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6]) #batch:批次
dim(expr)
#[1] 75376   170
anno=sample
dim(anno)
#234 6
anno=sample
dim(anno)
anno=anno[1:123,]
dim(anno)

expr1=t(expr)
expr1=expr1[1:123,]
dim(expr1)
expr1=t(expr1)
dim(expr1)
#[1]  75376 123
keep=rowSums(expr1>1) >= floor(0.05*ncol(expr1))
table(keep)
#FALSE  TRUE 
#18613 56763
expr <- expr1[keep,]
dim(expr)
expr1=t(expr)

#未去批次
library(ggplot2)
df_pca <- prcomp(log2(expr1+1)) #计算主成分
dim(expr1)
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)
head(df_pcs,3)
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+ geom_point()
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+
  geom_point(size=3)+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#D1B6E1","#fab1ce","#cbe86b","#4ea1d3"))


#去批次
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
dim(ad.batch)
expr1=expr1
class(expr1)
ad.clr=expr1

##removeBatchEffect
library(vegan)
library(limma)
library(mixOmics)
#removeBatchEffect
ad.mod <- model.matrix( ~ anno$batch)
dim(ad.clr)
ad.rBE <- t(removeBatchEffect(t(ad.clr), batch = anno$batch,design = ad.mod))

#PCA
library(gridExtra)
library("ggpubr")
library(ggplot2)
df_pca <- prcomp(ad.rBE) #计算主成分
dim(ad.rBE)
dim(expr1)
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)  

head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+ geom_point()
head(df_pcs,3)
plot(df_pca$x[,1], df_pca$x[,2])
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+ geom_point()
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+
  geom_point(size=3)+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#D1B6E1","#fab1ce","#cbe86b","#4ea1d3"))

#combat
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))
df_pca <- prcomp(ad.ComBat) #计算主成分
dim(ad.ComBat)
dim(expr1)
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=batch))+
  geom_point(size=3)+ 
  xlab(percentage[1]) +
  ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#D1B6E1","#fab1ce","#cbe86b","#4ea1d3"))




#######根（去批次）热图
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
library(xlsx)
sample <- read.xlsx("C:/Users/hr345/Desktop/all/各类sample/时间.xlsx",1,encoding="UTF-8")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6]) #batch:批次
dim(expr)
#[1] 75376   170

anno=sample
dim(anno)
anno=anno[1:123,]
dim(anno)
expr1=t(expr)
expr1=expr1[1:123,]
dim(expr1)
expr1=t(expr1)
dim(expr1)
#[1]  75376 123
keep=rowSums(expr1>1) >= floor(0.05*ncol(expr1))
table(keep)
#FALSE  TRUE 
#18613 56763
expr <- expr1[keep,]
dim(expr)
expr1=t(expr)


anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
anno$batch<- factor(anno$batch, levels=c("1","8","9","11"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("1","8","9","11"), ordered=TRUE)
dim(ad.batch)
ad.trt = factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
ad.factors.df <- data.frame(time = anno$time, batch = anno$batch)
expr1=expr1
class(expr1)
ad.clr=expr1
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))
expr2=ad.ComBat[,20000:25000]

expr1_unique=expr2
row.names(expr1_unique)=anno$time
sample_cor <- cor(t(expr1_unique))
sample_cor1 <- round(sample_cor, digits = 2)
#画图
library(pheatmap)
pheatmap(sample_cor1, display_numbers = F,fontsize = 10, angle_col = 0,
         cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,
		 cellheight = 3,cellweight = 8,fontsize_row = 5,fontsize_col = 5,
		 color = colorRampPalette(c("blue","white","red"))(100))



sample_dist <- dist(expr1_unique)
sample_hc <- hclust(sample_dist)
plot(sample_hc)

########筛选表达量
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
expr2=expr1[-171:-201,]
expr3=expr2[-48:-104,]
expr4=expr3[-47,]
expr5=expr4[-46,]
expr6=expr5[-32:-39,]
expr7=expr6[-25,]
expr8=expr7[-14:-22,]
expr9=expr8[-3:-4,]
expr10=expr9[-24,]
expr11=expr10[-18,]
expr12=expr11[-11,]
expr13=expr12[-3,]
expr=t(expr13)
expr=expr[-1,]
row.names(expr)= gsub(".[0-9].v2.1","",row.names(expr))
keep=rowSums(expr>1) >= floor(0.05*ncol(expr))
table(keep)
expr.all <- expr[keep,]
dim(expr.all)
known=fread("C:/Users/hr345/Desktop/all/1.txt",header=T,sep="\t")
cc=intersect(known$Geneid,row.names(expr.all))


##筛选TF
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/WGCNA/抗感/抗感sample.txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6])
dim(expr)
#[1] 75376   234
anno=sample
dim(anno)
#234 6
#####根做聚类   去掉1个样本
anno=sample[1:123,]
dim(anno)
expr1=t(expr)
expr1=expr1[1:123,]
expr1=t(expr1)
dim(expr1)
#[1]  75376 123
keep=rowSums(expr1>1) >= floor(0.05*ncol(expr1))
table(keep)
#FALSE  TRUE 
#18613 56763
expr <- expr1[keep,]
dim(expr)
#[1] 56763    123
row.names(expr)= gsub("[.][0-9].v2.1","",row.names(expr))
expr=t(expr)
anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
anno$concentration<- factor(anno$concentration, levels=c("0","150","200","salt-alkali","sk","mock"), ordered=TRUE)
anno$batch<- factor(anno$batch, levels=c("1","8","9","11"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("1","8","9","11"), ordered=TRUE)
dim(ad.batch)
ad.factors.df <- data.frame(time = anno$time, batch = anno$batch, concentration = anno$concentration)
expr1=expr
class(expr1)
ad.clr=expr1
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))
expr=t(ad.ComBat)
TF=fread("C:/Users/hr345/Desktop/all/TF/5048TF.txt",header=T,sep="\t")
write.csv(expr,file="expr.csv")#手动再第一行加上SRR，然后转成新建文本文档.txt
expr2=fread("C:/Users/hr345/Desktop/新建文本文档.txt",header=T,sep="\t")
aa=merge(TF,expr2,by="SRR",sort=F)
write.csv(aa,file="TF.csv")
library(dplyr)
specialcxy<-anti_join(expr2,aa,by='SRR')
write.csv(specialcxy,file="specialcxy.csv")#手动合并两个文件
sample2=t(anno)
write.csv(sample2,file="sample2.csv")#导出后将第一行SRR号提取出来粘贴到前面的文件中


