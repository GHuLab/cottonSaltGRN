########################################Leaf（Time+Variety）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample =fread("C:/Users/hr345/Desktop/all/sample/sample(final).txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6])
dim(expr)
#[1] 75376   170
anno=sample
dim(anno)
#234 6


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
df_pca <- prcomp(log2(expr1+1)) 
dim(expr1)
# 162 75376
anno$time<- factor(anno$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession)  
head(df_pcs,3) 
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") 
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

									   
###remove batch effect
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
ad.batch =  factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
dim(ad.batch)

ad.factors.df <- data.frame(batch = anno$batch)
expr1=expr1
class(expr1) 
ad.clr=expr1

#removeBatchEffect
#install.packages("vegan")
library(vegan)
library(limma)
library(mixOmics)
ad.mod <- model.matrix( ~ anno$batch)
dim(ad.clr)
ad.rBE <- t(removeBatchEffect(t(ad.clr), batch = anno$batch,design = ad.mod))
#PCA
library(gridExtra)
library("ggpubr")
library(ggplot2)
df_pca <- prcomp(ad.rBE) 
dim(ad.rBE)
#104 75376
dim(expr1)
#104 75376
anno$time<- factor(anno$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession)  
head(df_pcs,3) 
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf")
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


#####combat
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("sva")
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))
#PCA
df_pca <- prcomp(ad.ComBat)
dim(ad.ComBat)
#104 75376
dim(expr1)
#104 75376
anno$time<- factor(anno$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
anno$Accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,Accession=anno$Accession)  
head(df_pcs,3)
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf")
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
  scale_color_manual(values=c("#bbdffb","#90cbf9","#64b7f6","#41a7f5","#1e97f3","#1a8ae5","#1477d2","#1065c0","#0747a1","#09347A"))


###########Leaf（batch）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/sample/sample(final).txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6])
dim(expr)
#[1] 75376   170
anno=sample
dim(anno)
#234 6
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
##Do not temove batch effect
library(ggplot2)
df_pca <- prcomp(log2(expr1+1))
dim(expr1)
# 162 75376
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)  
head(df_pcs,3) 
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf")
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


#remove batch effect
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
df_pca <- prcomp(ad.rBE)
dim(ad.rBE)
#104 75376
dim(expr1)
#104 75376
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)  
head(df_pcs,3) 
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") 
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
df_pca <- prcomp(ad.ComBat)
dim(ad.ComBat)
#104 75376
dim(expr1)
#104 75376
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA722118","PRJNA253112","PRJNA248163","PRJNA490626"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)  

head(df_pcs,3)
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf")
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




#######Root（Time+Treatment）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/sample/sample(final).txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6]) 
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
df_pca <- prcomp(log2(expr1+1)) 
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
	

#remove batch effect
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
df_pca <- prcomp(ad.rBE)
dim(ad.rBE)
dim(expr1)
anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$concentration<- factor(anno$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,concentration=anno$concentration)  

head(df_pcs,3) 
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf")
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

#combat
library(sva)
ad.ComBat <- t(ComBat(t(ad.clr), batch = anno$batch))

#PCA
df_pca <- prcomp(ad.ComBat)
dim(ad.ComBat)
dim(expr1)
anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
anno$concentration<- factor(anno$concentration, levels=c("mock","150","200","salt-alkali","sk"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,concentration=anno$concentration)
head(df_pcs,3) 
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf")
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
     scale_color_manual(values=c("#ddead1","#c7ddb5","#b3cf99","#a3c585","#95bb72","#87ab69","#75975e","#658354","#4b6043"))+ scale_shape_manual(values=c(15,14,17,18,16))

#######Root（batch）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/E各类sample/sample(final).txt",header=T,sep="\t")
colnames(expr1)[1] <- 'SRR'
aa=merge(sample,expr1,by="SRR",sort=F)
dim(aa)
#[1]   234 75382
rownames(aa)=aa$SRR
expr=t(aa[,-1:-6])
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

#Do not remove batch effect
library(ggplot2)
df_pca <- prcomp(log2(expr1+1))
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


#remove batch effect
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
df_pca <- prcomp(ad.rBE)
dim(ad.rBE)
dim(expr1)
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)  

head(df_pcs,3) 
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf")
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
df_pca <- prcomp(ad.ComBat)
dim(ad.ComBat)
dim(expr1)
anno$batch<- factor(anno$batch, levels=c("PRJNA531727","PRJNA485838","PRJNA532694","PRJNA919499"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,batch = anno$batch)
head(df_pcs,3) 
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") 
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




#######Root's pheatmeat
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
expr=t(aa[,-1:-6])
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

library(pheatmap)
pheatmap(sample_cor1, display_numbers = F,fontsize = 10, angle_col = 0,
         cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,
		 cellheight = 3,cellweight = 8,fontsize_row = 5,fontsize_col = 5,
		 color = colorRampPalette(c("blue","white","red"))(100))



sample_dist <- dist(expr1_unique)
sample_hc <- hclust(sample_dist)
plot(sample_hc)

