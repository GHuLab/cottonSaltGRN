##叶片（按批次分）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
library(xlsx)
sample <- read.xlsx("C:/Users/hr345/Desktop/all/各类sample/sample（去批次）2.xlsx",1,encoding="UTF-8")
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

library(ggplot2)
library(umap)
iris.umap = umap::umap(expr1)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$batch)
colnames(pumap) = c("UMAP_1","UMAP_2","batch")
rownames(pumap)=anno$SRR
head(pumap)
pumap$batch<- factor(pumap$batch, levels=c("PRJNA531727","PRJNA601953","PRJNA623201","PRJNA253112","PRJNA722118","PRJNA248163","PRJNA490626"), ordered=TRUE)
ggplot(pumap,aes(UMAP_1,UMAP_2,color=batch))+
       geom_point(size=3)+ 
       theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
                         legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
       theme(axis.title.y = element_text(size = 15))+ 
       scale_color_manual(values=c("#a3daff","#0080ff","#4ea1d3","#ede574","#fc913a","#ff4e50","#cbe86b"))

##叶片（按种类分）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
library(xlsx)
sample <- read.xlsx("C:/Users/hr345/Desktop/all/各类sample/sample（去批次）2.xlsx",1,encoding="UTF-8")
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

library(ggplot2)
library(umap)
iris.umap = umap::umap(expr1)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$Accession)
colnames(pumap) = c("UMAP_1","UMAP_2","Accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$Accession<- factor(pumap$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
ggplot(pumap,aes(UMAP_1,UMAP_2,color=Accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#fc913a","#cbe86b"))


##叶片（按处理分）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
library(xlsx)
sample <- read.xlsx("C:/Users/hr345/Desktop/all/各类sample/sample（去批次）2.xlsx",1,encoding="UTF-8")
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

library(ggplot2)
library(umap)
iris.umap = umap::umap(expr1)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$concentration)
colnames(pumap) = c("UMAP_1","UMAP_2","concentration")
rownames(pumap)=anno$SRR
head(pumap)
pumap$concentration<- factor(pumap$concentration, levels=c("0","salt-alkali","50","400","200","100","un"), ordered=TRUE)
ggplot(pumap,aes(UMAP_1,UMAP_2,color=concentration))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#0080ff","#4ea1d3","#ede574","#fc913a","#ff4e50","#cbe86b"))


##叶片（按时间分）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
library(xlsx)
sample <- read.xlsx("C:/Users/hr345/Desktop/all/各类sample/sample（去批次）2.xlsx",1,encoding="UTF-8")
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

library(ggplot2)
library(umap)
iris.umap = umap::umap(expr1)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$time)
colnames(pumap) = c("UMAP_1","UMAP_2","time")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h", "1h","3h","4h","6h","12h","24h","48h","2w","mock"), ordered=TRUE)
ggplot(pumap,aes(UMAP_1,UMAP_2,color=time))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#a3daff","#1ec0ff","#0080ff","#4ea1d3","#ede574","#f8ca00","#fc913a","#ff4e50","#cbe86b","#3b8686"))

##根（按批次分）
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
library(xlsx)
sample <- read.xlsx("C:/Users/hr345/Desktop/all/各类sample/sample（去批次）2.xlsx",1,encoding="UTF-8")
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

#UMAP
library(ggplot2)
library(umap)
iris.umap = umap::umap(expr1)
head(iris.umap$layout)
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$batch)
colnames(pumap) = c("UMAP_1","UMAP_2","batch")
rownames(pumap)=anno$SRR
head(pumap)
pumap$batch<- factor(pumap$batch, levels=c("PRJNA531727", "PRJNA485838", "PRJNA532694","PRJNA919499"), ordered=TRUE)


ggplot(pumap,aes(UMAP_1,UMAP_2,color=batch))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#D1B6E1","#EE7785","#8EC0E4","#cff09e"))
