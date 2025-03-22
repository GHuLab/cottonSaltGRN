rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=read.table("C:/Users/hr345/Desktop/all/allsample.txt",header=T,sep="\t")
colnames(expr1)[1] <- 'sample'#将sample的第一列变成expr1的第一列
aa=merge(sample,expr1,by="sample",sort=F)
dim(aa)
#[1]   234 75381
rownames(aa)=aa$SRR
expr=t(aa[,-1:-5])
dim(expr)
#[1] 75376   234

######耐盐的基因###
keep=rowSums(expr>1) >= floor(0.5*ncol(expr))
table(keep)
#keep
#FALSE  TRUE 
#30458 44918 
expr <- expr[keep,]
dim(expr)
#[1] 44918   234

anno=sample
dim(anno)
#234 5


#########根+叶

anno=anno[-92,]
anno=anno[-162:-167,]
dim(anno)
expr1=t(expr)
expr1=expr1[-92,]
expr1=expr1[-162:-167,]
dim(expr1)
expr1_unique<-unique(expr1)
dim(expr1_unique)
class(expr1_unique)

#树状图
dists <- dist(expr1,method = "euclidean")
hc <- hclust(dists, method = "ave")
dend1 <- as.dendrogram(hc)
plot(dend1, type = "rectangle", 
     ylab="Height",
     main="Cluster Dendrogram")

#tSNE
library(Rtsne)
tsne_out = Rtsne(expr1_unique,perplexity=20)
str(tsne_out)
pdat = data.frame(tsne_out$Y)
rownames(pdat)=anno$SRR
pdat=cbind(pdat,anno$tissue)
colnames(pdat) = c("tSNE_1","tSNE_2","tissue")
head(pdat)
pdat$tissue<- factor(pdat$tissue, levels=c("R","L"), ordered=TRUE)
library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=tissue,shape=tissue))+
  geom_point(size=3,alpha=1/2)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#6a60a9","#f199bc"))
  
#PCA
library(ggplot2)
df_pca <- prcomp(log2(expr1+1)) #计算主成分
dim(expr1)
anno$tissue<- factor(anno$tissue, levels=c("R","L"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,tissue = anno$tissue)
head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=tissue,shape=tissue))+ geom_point() 
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=tissue,shape=tissue))+
  geom_point(size=3, alpha=1/2)+ 
  # xlab(percentage[1]) +
  # ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#6a60a9","#f199bc")) 
  
#UMAP
library(umap)
#expr1=t(expr)
dim(expr1)
expr1_unique<-unique(expr1)
dim(expr1_unique)
# 使用umap函数进行UMAP降维分析
iris.umap = umap::umap(expr1_unique)
# 查看降维后的结果
head(iris.umap$layout)
##           [,1]     [,2]
#SRR10426360 -3.447716 1.967124
#SRR10426373 -3.613967 2.087851
#SRR12077547 -4.925000 3.021000
#SRR12077573 -4.894924 2.986124
# 使用plot函数可视化UMAP的结果]
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$tissue)
colnames(pumap) = c("UMAP_1","UMAP_2","tissue")
rownames(pumap)=anno$SRR
head(pumap)
pumap$tissue<- factor(pumap$tissue, levels=c("R","L"), ordered=TRUE)
ggplot(pumap,aes(UMAP_1,UMAP_2,color=tissue,shape=tissue))+
  geom_point(size=3, alpha=1/2)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#6a60a9","#f199bc"))


#####根做聚类   去掉1个样本

anno=anno[-88:-198,]
dim(anno)
#123 5
#T-SNE  tissue
library(Rtsne)
expr1=t(expr)
expr1=expr1[-88:-198,]
dim(expr1)
#[1]   123 44918
#T-SNE  tissue
#第三次
#87个样本  浓度
expr1_unique<-unique(expr1)
dim(expr1_unique)
#[1]    123 44918
class(expr1_unique)

tsne_out = Rtsne(expr1_unique,perplexity=20)
str(tsne_out)
pdat = data.frame(tsne_out$Y)
rownames(pdat)=anno$SRR
pdat=cbind(pdat,anno$concentration,anno$Accession)

colnames(pdat) = c("tSNE_1","tSNE_2","concentration","accession")
head(pdat)
pdat$concentration<- factor(pdat$concentration, levels=c("0","150","200","salt-alkali","s0","sk","mock"), ordered=TRUE)
pdat$accession<- factor(pdat$accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=concentration,shape=accession))+
  geom_point(size=3,alpha=1/2)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#f9a11b","#d62a9d","#6a60a9","#f199bc"))
#PCA
###PCA analysis
library(ggplot2)
df_pca <- prcomp(log2(expr1+1)) #计算主成分
dim(expr1)
# 87 45140
anno$concentration<- factor(anno$concentration, levels=c("0","150","200","salt-alkali","s0","sk","mock"), ordered=TRUE)
anno$accession<- factor(anno$Accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,concentration = anno$concentration,accession=anno$Accession)  

head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=concentration,shape=accession))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,color=concentration,shape=accession))+
  geom_point(size=3, alpha=1/2)+ 
  # xlab(percentage[1]) +
  # ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#f9a11b","#d62a9d","#6a60a9","#f199bc"))
#700 500

#UMAP
library(umap)
#expr1=t(expr)
dim(expr1)
expr1_unique<-unique(expr1)
dim(expr1_unique)
# 使用umap函数进行UMAP降维分析
iris.umap = umap::umap(expr1_unique)
# 查看降维后的结果
head(iris.umap$layout)
##           [,1]     [,2]
#SRR10426360 -3.447716 1.967124
#SRR10426373 -3.613967 2.087851
#SRR12077547 -4.925000 3.021000
#SRR12077573 -4.894924 2.986124
# 使用plot函数可视化UMAP的结果]
pumap = data.frame(iris.umap$layout)

pumap=cbind(pumap,anno$concentration,anno$accession)
colnames(pumap) = c("UMAP_1","UMAP_2","concentration","accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$concentration<- factor(pumap$concentration, levels=c("0","150","200","salt-alkali","s0","sk","mock"), ordered=TRUE)
pumap$accession<- factor(pumap$accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=concentration,shape=accession))+
  geom_point(size=3, alpha=1/2)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ scale_color_manual(values=c("#f9320c","#60c5ba","#03a6ff","#f9a11b","#d62a9d","#6a60a9","#f199bc"))




#时间
tsne_out = Rtsne(expr1_unique,perplexity=30)
pdat = data.frame(tsne_out$Y)
pdat=cbind(pdat,anno$time,anno$accession)
colnames(pdat) = c("tSNE_1","tSNE_2","time","accession")
head(pdat)
pdat$time<- factor(pdat$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h",
                                       "72h","s0h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
pdat$accession<- factor(pdat$accession, levels=c("T","S","TM-1"), ordered=TRUE)

library(ggplot2)
ggplot(pdat,aes(x=tSNE_1,y=tSNE_2,color=time,shape=accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#77919d","#5e5e5f","#88dba3","#1f4e5f","#cbe86b","#4f953b",
                                       "#a3daff","#1ec0ff","#0080ff","#4ea1d3","#eb9f9f","#bf209f","#6a60a9","#5c196b"))
									   
#99FFFF
#PCA
library(ggplot2)
#df_pca <- prcomp(log2(expr1+1)) #计算主成分
dim(expr1)
# 162 75376
anno$time<- factor(anno$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h",
                                       "72h","s0h","s3h","s12h","s48h","sk1h","sk6h","sk12h","sk24h"), ordered=TRUE)
anno$accession<- factor(anno$accession, levels=c("T","S","TM-1"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,time = anno$time,accession=anno$accession)  

head(df_pcs,3)  #查看主成分结果
plot(df_pca$x[,1], df_pca$x[,2])
#pdf("pca.log2norm413.pdf") #可以单独画画
ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=accession))+ geom_point()

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

ggplot(df_pcs,aes(x=PC1,y=PC2,color=time,shape=accession))+
  geom_point(size=3)+ 
  # xlab(percentage[1]) +
  # ylab(percentage[2])+
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#77919d","#5e5e5f","#88dba3","#60c5ba","#cbe86b","#4f953b",
                                       "#a3daff","#1ec0ff","#0080ff","#4ea1d3","#eb9f9f","#bf209f","#6a60a9","#5c196b"))
#700 500

#UMAP
library(umap)
#expr1=t(expr)
dim(expr1)
expr1_unique<-unique(expr1)
dim(expr1_unique)
# 使用umap函数进行UMAP降维分析
#iris.umap = umap::umap(expr1_unique)
# 查看降维后的结果
head(iris.umap$layout)
##           [,1]     [,2]
#SRR10426360 -3.447716 1.967124
#SRR10426373 -3.613967 2.087851
#SRR12077547 -4.925000 3.021000
#SRR12077573 -4.894924 2.986124
# 使用plot函数可视化UMAP的结果]
pumap = data.frame(iris.umap$layout)
pumap=cbind(pumap,anno$time,anno$accession)
colnames(pumap) = c("UMAP_1","UMAP_2","time","accession")
rownames(pumap)=anno$SRR
head(pumap)
pumap$time<- factor(pumap$time, levels=c("0h", "15m", "1h","3h","6h","12h","24h","48h",
                                         "72h","s0h","s3h","s12h","s48h"), ordered=TRUE)
pumap$accession<- factor(pumap$accession, levels=c("T","S","TM-1"), ordered=TRUE)

ggplot(pumap,aes(UMAP_1,UMAP_2,color=time,shape=accession))+
  geom_point(size=3)+ 
  theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
        legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15))+ 
  scale_color_manual(values=c("#fbd14b","#f26d5b","#c03546","#77919d","#5e5e5f","#88dba3","#60c5ba","#cbe86b","#4f953b",
                                       "#a3daff","#1ec0ff","#0080ff","#4ea1d3","#eb9f9f","#bf209f","#6a60a9","#5c196b"))
									   

#相关性分析
#计算距离

expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")

colnames(expr1)[1] <- 'sample'
sample=read.table("C:/Users/hr345/Desktop/all/allsample.txt",header=T,sep="\t")
aa=merge(sample,expr1,by="sample",sort=F)

dim(aa)
#167 75381
rownames(aa)=aa$SRR
expr=t(aa[,-1:-5])
dim(expr)
#[1] 75376   167
anno=sample
dim(anno)
#167 5

#####根做聚类   去掉1个样本
anno=sample[-88:-198,]

dim(anno)
#87 5
#T-SNE  tissue
library(Rtsne)
expr1=t(expr)
expr1=expr1[-88:-198,]


dim(expr1)

expr1_unique=expr1

row.names(expr1_unique)=anno$time
sample_cor <- cor(t(expr1_unique))
sample_cor1 <- round(sample_cor, digits = 2)
#画图
library(pheatmap)
pheatmap(sample_cor1, display_numbers = F,fontsize = 10, angle_col = 90,
         cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = T,
		 cellheight = 5,cellweight = 5,fontsize_row = 10,fontsize_col = 10)



sample_dist <- dist(expr1_unique)
sample_hc <- hclust(sample_dist)
plot(sample_hc)





########WGCNA

###########根
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=read.table("C:/Users/hr345/Desktop/all/allsample.txt",header=T,sep="\t")
colnames(expr1)[1] <- 'sample'
aa=merge(sample,expr1,by="sample",sort=F)
dim(aa)
rownames(aa)=aa$SRR
expr=t(aa[,-1:-5])
dim(expr)
#[1] 75376   234
anno=sample
dim(anno)
#234 5
#####根做聚类   去掉1个样本
anno=sample[-88:-198,]

dim(anno)
#123 5
#T-SNE  tissue
library(Rtsne)
expr1=t(expr)
expr1=expr1[-88:-198,]



expr1=t(expr1)
dim(expr1)
#[1]  75376 123
keep=rowSums(expr1>1) >= floor(0.5*ncol(expr1))
table(keep)
#FALSE  TRUE 
#28390 46986
expr <- expr1[keep,]
dim(expr)
#[1] 46986    87
row.names(expr)= gsub("[.][0-9].v2.1","",row.names(expr))
TF=read.table('C:/Users/hr345/Desktop/all/tf.txt',header=F,sep="\t")
#####一列基因ID
cc1=intersect(TF$V1,row.names(expr))
#50% >0   TF  4047/5036   gene    59710/75376

#setdiff(TF$V1,row.names(expr))
srr12=anno
coldata1<-data.frame(srr=srr12$sample,sample = srr12$sample)
coldata=coldata1
head(coldata)
#srr      sample
#1 SRR10518033 SRR10518033
#2 SRR10518014 SRR10518014
dim(coldata)
#123 2
library(RColorBrewer)
library(WGCNA)
#install.packages("backports")
#library(backports)

library(flashClust)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

multiExpr = vector(mode = "list", length = 1)
multiExpr[[1]] = list(data = t(expr));


checkSets(multiExpr)
# $nSets [1] 1
# $nGenes [1]46986
# $nSamples [1] 123
# $structureOK [1] TRUE
nSets=checkSets(multiExpr)$nSets

filterMultiExpr<-function(multiExpr,nSets)
{ gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
print(gsg$allOK)
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
  for (set in 1:nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  } }
print(checkSets(multiExpr))
return(multiExpr)
}
multiExpr<-filterMultiExpr(multiExpr,nSets)
# $nSets [1] 1
# $nGenes [1] 46986
# $nSamples [1] 123
# $structureOK [1] TRUE
shortLabels="All"
save(expr, coldata, multiExpr,  file="C:/Users/熊显鹏/桌面/salt/Rtemp59712.rdata")

## Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
# types<-c("unsigned", "signed", "signed hybrid")
powers = c(c(1:10), seq(from = 12, to=40, by=2))
for(type in c("signed", "signed hybrid")){
  for(corU in c("cor","bicor")){
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
      powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, corFnc = get(corU), networkType = type, blockSize=10000)[[2]])      }
    collectGarbage()
    # Plot the results:
    colors=brewer.pal(nSets,"Set1")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
      for (col in 1:length(plotCols))
      {
        ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
      }
    }
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    # sizeGrWindow(8, 6)
    pdf(paste0("wgcna.choosePower.",gsub(" ","",type),".",corU,".pdf") )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
      if (set==1)
      {
        plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
        addGrid()
      }
      if (col==1)
      {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
             labels=powers,cex=cex1,col=colors[set]);
      } else
        text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
      if (col==1)
      {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
      } else
        legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
    dev.off()
    assign(paste0("powerTables.",gsub(" ","",type),".",corU) , powerTables)
  }
}
# examine plot and also print out power choice, plus default 12

powerTables=list()
for(type in c("unsigned", "signed", "signed hybrid"))
{
  for(corU in c("cor","bicor")){
    p=paste0("powerTables.",gsub(" ","",type),".",corU)
    print(p)
    print(lapply(get(p), function(x){ x<-x$data; x$Power[x$SFT.R.sq>0.8 & x$slope<0][1]} ))
    powerTables[[paste0(gsub(" ","",type),".",corU)]] <-get(p)[[1]]
  }
}

powers=c(12,16,24)
save(expr, multiExpr, powerTables, powers, file = "C:/Users/熊显鹏/桌面/salt/wgcna.TPM.59712.prep.Rdata")


for(corM in c("bicor","pearson")){
  for(power in powers){
    # construct A2D5 network
    cgn =  blockwiseModules(
      # Input data
      multiExpr[[1]]$data,
      # Data checking options
      checkMissingData = TRUE,
      # Options for splitting data into blocks
      maxBlockSize =  30000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
      #randomSeed = 12345,
      # Network construction arguments: correlation options, use bicor instead of default pearson
      corType = corM,
      # Adjacency and topology overlap function options
      power = power, networkType = "signed", TOMType = "signed",
      # load previous TOMs
      saveTOMs = FALSE,
      # Basic tree cut options
      deepSplit = 2,  #default, known to reasonable
      minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
      pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
      # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
      mergeCutHeight = 0.25,
      # others
      reassignThreshold = 0,
      numericLabels = TRUE,
      verbose = 3)
    assign(paste0("cgnP",power,substring(corM,1,1)),cgn)
    save(list=grep("cgnP",ls(),value=TRUE), file = "wgcna.cgn1104.Rdata")
  }
}




# 34672 genes    直接用这个跑后面的
library(RColorBrewer)
library(flashClust)
library(WGCNA)
load("C:/Users/熊显鹏/桌面/salt/root87WGCNA/34672/Rtemp34672.rdata")
load("C:/Users/熊显鹏/桌面/salt/root87WGCNA/34672/wgcna.34672cgn.Rdata")
load("C:/Users/熊显鹏/桌面/salt/root87WGCNA/34672/wgcna.TPM.34672.prep.Rdata")

ls()
i="cgnP30p"
net = get(i)
datExpr=as.data.frame(t(expr))
dim(datExpr)
# [1]413 25441
moduleLabels = net$colors

moduleColors = labels2colors(net$colors)
table(moduleColors)

#moduleColors
#black         blue        brown         cyan        green  greenyellow         grey    lightcyan 
#624         2565         1499          187         1160          291        16059          128 
#magenta midnightblue         pink       purple          red       salmon          tan    turquoise 
#453          156          504          328         1059          196          277         7742 
#yellow 
#1444 
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
dim(MET)
#87 17
pdf(file = "Eigengene adjacency heatmap.pdf", width = 12, height = 9)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


#第二种
datTraits=read.table("C:/E盘/博士后第一个项目/最新的1019/20211027TPM/1110/Datatraits1.txt",header=T,sep="\t",row.names = 1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
moduleTraitCor1 = cor(MEs, datTraits, use = "p")
moduleTraitPvalue1 = corPvalueStudent(moduleTraitCor1, nSamples)
textMatrix1 = paste(signif(moduleTraitCor1, 2), "\n(", signif(moduleTraitPvalue1, 1), ")", sep = "")
dim(textMatrix1)
dim(moduleTraitCor1)
par(mar = c(5, 10, 3, 3))
#datTraits1=data.frame(datTraits[,-2])
#names(datTraits1)="DPA"
labeledHeatmap(Matrix = moduleTraitCor1, xLabels = names(datTraits), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix1, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), 
               main = paste("Module-Fiber development period relationships"))



# WGCNA hub gene heatmap
MEs1=t(MEs)
c=data.frame(cbind(apply(MEs1[,5:7],1,mean),
                   apply(MEs1[,11:13],1,mean),
                   apply(MEs1[,20:22],1,mean),
                   apply(MEs1[,28:30],1,mean),
                   apply(MEs1[,34:36],1,mean),
                   apply(MEs1[,46:48],1,mean),
                   apply(MEs1[,55:57],1,mean),
                   apply(MEs1[,67:69],1,mean),
                   apply(MEs1[,8:10],1,mean),
                   apply(MEs1[,14:16],1,mean),
                   apply(MEs1[,17:19],1,mean),
                   apply(MEs1[,31:33],1,mean),
                   apply(MEs1[,37:39],1,mean),
                   apply(MEs1[,49:51],1,mean),
                   apply(MEs1[,52:54],1,mean),
                   apply(MEs1[,64:66],1,mean),
                   MEs1[,1],
                   apply(MEs1[,23:24],1,mean),
                   apply(MEs1[,28:29],1,mean),
                   apply(MEs1[,58:59],1,mean),
                   apply(MEs1[,2:4],1,mean),
                   apply(MEs1[,25:27],1,mean),
                   apply(MEs1[,42:45],1,mean),
                   apply(MEs1[,60:63],1,mean),
                   apply(MEs1[,70:72],1,mean),
                   apply(MEs1[,73:75],1,mean),
                   apply(MEs1[,76:78],1,mean),
                   apply(MEs1[,79:81],1,mean),
                   apply(MEs1[,82:84],1,mean),
                   apply(MEs1[,85:87],1,mean)))

class(c)
#[1] "data.frame"
#修改行名
names(c) <- c("a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o1",
              "a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1")

#c$module=names(MEs)
c$module=c("black(624)","blue(2565)","brown(1499)","cyan(187)","green(624)",
           "greenyellow(2565)","grey(16059)",
           "lightcyan(128)","magenta(453)", "midnightblue(156)","pink(504)",
           "purple(328)","red(1069)","salmon(196)","tan(277)",
           "turquoise(7742)","yellow(1444)")



c <- c[, c("module","a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o1",
           "a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1")]

anno=read.table("C:/Users/熊显鹏/桌面/salt/root87/34672/87.Sample.information.2.txt",header=T,sep="\t")
Period=anno[,-1]
#names(Period)="time"
head(Period)
#Period  DPA
#1 Elongation 5dpa
class(Period)
dim(Period)
#anno11=data.frame(Period) 非数值需要转置


#per=read.table("C:/E盘/博士后第一个项目/最新的1019/20211027TPM/1110/0420/热图/MEs.4.txt",header=T, sep="\t")
per=c
dat2 <- per[,-1]*10
rownames(Period) = rownames(t(dat2))
row.names(dat2)=per$module
#dat2=t(dat2)
library(pheatmap)
#Period$time <- factor(Period$time, levels=c("S0", "S15min", "S1", "S3", "S6", "S12", "S24", "S48", 
#                                            "T0", "T15min", "T1", "T3", "T6", "T12", "T24", "T48", 
#                                           "saS0","saS3","saS12","saS48",
#                                          "saT0","saT3","saT12","saT48",
#                                       "TM0h", "TM6h","TM12h","TM24h",
#                                      "TM48h","TM72h"), ordered=TRUE)

#ann_colors1 = list(DPA = c("5dpa"="#008EFF", "7dpa"="#00BAFF", "8dpa"="#00E5FF", "10dpa"="#00FF4D",
#                      "12dpa"="#6BFF00", "13dpa"="#A8FF00", "15dpa"="#FFFF80", 
#                     "18dpa"="#FFFF80", "19dpa"="#FFFF00", "20dpa"="#FFE500", "24dpa"="#FFB300", 
#                    "25dpa"="#FF8000", "28dpa"="#FF4D00", "30dpa"="#FF0000"))

pheatmap(dat2[,1:16],annotation_col = Period,
         #annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,
         #fontface="italic",
         fontfamily= "Times New Roman",show_colnames = F)

pheatmap(dat2,annotation_col = Period,
         #annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,
         #fontface="italic",
         fontfamily= "Times New Roman",show_colnames = F,gaps_col = c(8, 16,20,24))


dat2=dat2[-7,-25:-30]
dat2=dat2[,-15]
dat2=dat2[,-13]
dat2=dat2[,-10:-11]
dat2=dat2[,-7]
dat2=dat2[,-5]
dat2=dat2[,-2:-3]
Period$Accession<- factor(Period$Accession, levels=c("susceptible(salt)", "tolerant(salt)",
                                                     "susceptible(salt-alkali)","tolerant(salt-alkali)"), ordered=TRUE)
Period$time<- factor(Period$time, levels=c("0h", "3h", "12h","48h"), ordered=TRUE)

Period$time <- factor(Period$time, levels=c("S0",  "S3", "S12", "S48", 
                                            "T0", "T3", "T12", "T48", 
                                            "saS0","saS3","saS12","saS48",
                                            "saT0","saT3","saT12","saT48"), ordered=TRUE)

ann_colors1 = list(Accession = c("susceptible(salt)"="#7FCDBB", "tolerant(salt)"= "#41B6C4",
                                 "susceptible(salt-alkali)"="#1D91C0","tolerant(salt-alkali)"="#225EA8"),
                   time = c("0h"="#FFF7BC",  "3h"="#FEE391", "12h"= "#FEC44F", "48h"="#FE9929"))




pheatmap(dat2,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,
         #fontface="italic",
         fontfamily= "Times New Roman",show_colnames = F,
         gaps_col = c(4, 8,12,16))







library(WGCNA);
library(RColorBrewer)
library(flashClust);
load("C:/Users/熊显鹏/桌面/salt/root87/34672/Rtemp34672.rdata")
load("C:/Users/熊显鹏/桌面/salt/root87/34672/wgcna.34672cgn.Rdata")
load("C:/Users/熊显鹏/桌面/salt/root87/34672/wgcna.TPM.34672.prep.Rdata")

ls()
i="cgnP30p"
net = get(i)
datExpr=as.data.frame(t(expr))
dim(datExpr)
# [1]413 25441
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)

#moduleColors
#black         blue        brown         cyan        green  greenyellow         grey    lightcyan 
#624         2565         1499          187         1160          291        16059          128 
#magenta midnightblue         pink       purple          red       salmon          tan    turquoise 
#453          156          504          328         1059          196          277         7742 
#yellow 
#1444 
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
dim(MET)

# WGCNA hub gene heatmap
MEs1=t(MEs)
c=data.frame(cbind(apply(MEs1[,5:7],1,mean),
                   apply(MEs1[,11:13],1,mean),
                   apply(MEs1[,20:22],1,mean),
                   apply(MEs1[,28:30],1,mean),
                   apply(MEs1[,34:36],1,mean),
                   apply(MEs1[,46:48],1,mean),
                   apply(MEs1[,55:57],1,mean),
                   apply(MEs1[,67:69],1,mean),
                   apply(MEs1[,8:10],1,mean),
                   apply(MEs1[,14:16],1,mean),
                   apply(MEs1[,17:19],1,mean),
                   apply(MEs1[,31:33],1,mean),
                   apply(MEs1[,37:39],1,mean),
                   apply(MEs1[,49:51],1,mean),
                   apply(MEs1[,52:54],1,mean),
                   apply(MEs1[,64:66],1,mean),
                   MEs1[,1],
                   apply(MEs1[,23:24],1,mean),
                   apply(MEs1[,28:29],1,mean),
                   apply(MEs1[,58:59],1,mean),
                   apply(MEs1[,2:4],1,mean),
                   apply(MEs1[,25:27],1,mean),
                   apply(MEs1[,42:45],1,mean),
                   apply(MEs1[,60:63],1,mean),
                   apply(MEs1[,70:72],1,mean),
                   apply(MEs1[,73:75],1,mean),
                   apply(MEs1[,76:78],1,mean),
                   apply(MEs1[,79:81],1,mean),
                   apply(MEs1[,82:84],1,mean),
                   apply(MEs1[,85:87],1,mean)))

class(c)
#[1] "data.frame"
#修改行名
names(c) <- c("a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o1",
              "a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1")

#c$module=names(MEs)
c$module=c("black(624)","blue(2565)","brown(1499)","cyan(187)","green(624)",
           "greenyellow(2565)","grey(16059)",
           "lightcyan(128)","magenta(453)", "midnightblue(156)","pink(504)",
           "purple(328)","red(1069)","salmon(196)","tan(277)",
           "turquoise(7742)","yellow(1444)")



c <- c[, c("module","a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o1",
           "a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1")]

anno=read.table("C:/Users/熊显鹏/桌面/salt/root87/34672/87.Sample.information.2.txt",header=T,sep="\t")
Period=anno[,-1]
#names(Period)="time"
head(Period)
#Period  DPA
#1 Elongation 5dpa
class(Period)
dim(Period)
#anno11=data.frame(Period) 非数值需要转置


#per=read.table("C:/E盘/博士后第一个项目/最新的1019/20211027TPM/1110/0420/热图/MEs.4.txt",header=T, sep="\t")
per=c
dat2 <- per[,-1]*10
rownames(Period) = rownames(t(dat2))
row.names(dat2)=per$module
#dat2=t(dat2)
library(pheatmap)
#Period$time <- factor(Period$time, levels=c("S0", "S15min", "S1", "S3", "S6", "S12", "S24", "S48", 
#                                            "T0", "T15min", "T1", "T3", "T6", "T12", "T24", "T48", 
#                                           "saS0","saS3","saS12","saS48",
#                                          "saT0","saT3","saT12","saT48",
#                                       "TM0h", "TM6h","TM12h","TM24h",
#                                      "TM48h","TM72h"), ordered=TRUE)

#ann_colors1 = list(DPA = c("5dpa"="#008EFF", "7dpa"="#00BAFF", "8dpa"="#00E5FF", "10dpa"="#00FF4D",
#                      "12dpa"="#6BFF00", "13dpa"="#A8FF00", "15dpa"="#FFFF80", 
#                     "18dpa"="#FFFF80", "19dpa"="#FFFF00", "20dpa"="#FFE500", "24dpa"="#FFB300", 
#                    "25dpa"="#FF8000", "28dpa"="#FF4D00", "30dpa"="#FF0000"))

pheatmap(dat2[,1:16],annotation_col = Period,
         #annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,
         #fontface="italic",
         fontfamily= "Times New Roman",show_colnames = F)

pheatmap(dat2[,-25:-30],annotation_col = Period,
         #annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,
         #fontface="italic",
         fontfamily= "Times New Roman",show_colnames = F,gaps_col = c(8, 16,20,24))


dat2=dat2[-7,-25:-30]
#dat2=dat2[,-15]
#dat2=dat2[,-13]
#dat2=dat2[,-10:-11]
#dat2=dat2[,-7]
#dat2=dat2[,-5]
#dat2=dat2[,-2:-3]
Period$Accession<- factor(Period$Accession, levels=c("susceptible(salt)", "tolerant(salt)",
                                                     "susceptible(salt-alkali)","tolerant(salt-alkali)"), ordered=TRUE)
Period$time<- factor(Period$time, levels=c("0h", "15min","1h","3h", 
                                           "6h","12h","24h","48h"), ordered=TRUE)

ann_colors1 = list(Accession = c("susceptible(salt)"="#7FCDBB", "tolerant(salt)"= "#41B6C4",
                                 "susceptible(salt-alkali)"="#1D91C0","tolerant(salt-alkali)"="#225EA8"),
                   time = c("0h"="#FFFF80","15min"="#00FF4D", "1h"="#A8FF00" ,
                            "3h"="#00E5FF", "6h"="#008EFF", "12h"="#FFB300", 
                            "24h"="#FF8000","48h"="#FF4D00"))



pheatmap(dat2,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,
         #fontface="italic",
         fontfamily= "Times New Roman",show_colnames = F,
         gaps_col = c(8,16,20))











#pdf 导出
library(RColorBrewer)
display.brewer.all() 

brewer.pal(9,"YlOrBr")

display.brewer.pal(9,"YlGnBu")

