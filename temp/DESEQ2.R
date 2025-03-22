rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/counts/PRJNA532649/counts.txt",header=T,sep="\t")
OA=fread("C:/Users/hr345/Desktop/counts/OA.txt",header=T,sep="\t")
OD=fread("C:/Users/hr345/Desktop/counts/OD.txt",header=T,sep="\t")
colnames(expr1)[1] <- 'OA'
aa=merge(expr1,OA,by="OA",sort=F)
colnames(expr1)[1] <- 'OD'
bb=merge(expr1,OD,by="OD",sort=F)
write.csv(aa,file="C:/Users/hr345/Desktop/counts/PRJNA532649/OA.csv")
write.csv(bb,file="C:/Users/hr345/Desktop/counts/PRJNA532649/OD.csv")

####Deseq2
library(DESeq2)
#显示mycounts信息
rm(list=ls())
#设置工作目录
setwd("C:/Users/hr345/Desktop/all/counts")
exprSetTPM=read.table('C:/Users/hr345/Desktop/all/counts/PRJNA919499/count4.txt',header=T,sep="\t",row.names=1)
head(exprSetTPM)
#WT_rep1 WT_rep2 WT_rep3 GhWRKY16.RNAi_rep1 GhWRKY16.RNAi_rep2 GhWRKY16.RNAi_rep3 GhWRKY16.RNAi_rep5
#Gohir.A01G000101       1       6       6               1                4                4                 11
#Gohir.A01G000201      18      21      17               15               14               13                 54
###一般3个重复就行，这里有4个，所以我删除了1个
exprSetTPM1=exprSetTPM[,-4:-36]
#exprSetTPM2=exprSetTPM1[,-4:-36]
#exprSetTPM2=exprSetTPM1[,-4:-48]
exprSetTPM2=exprSetTPM1[,-4:-18]
exprSetTPM3=exprSetTPM2[,1:6]
head(exprSetTPM3)
#   WT_rep1 WT_rep2 WT_rep3 GhWRKY16.RNAi_rep2 GhWRKY16.RNAi_rep3 GhWRKY16.RNAi_rep5
#Gohir.A01G000101       1       6       6                  4                  4                 11
#Gohir.A01G000201      18      21      17                 14                 13                 54

##设置样品组别、重复数
#condition <- factor(c(rep("CK", 3), rep("CK1",3)),levels = c("CK","CK1"))
#condition <- factor(c(rep("SK", 3), rep("SK1",3)),levels = c("SK","SK1"))
#condition <- factor(c(rep("S", 3), rep("S1",3)),levels = c("S","S1"))
#condition <- factor(c(rep("H15", 3), rep("DH15",3)),levels = c("H15","DH15"))
condition <- factor(c(rep("ZM12", 3), rep("DZM12",3)),levels = c("ZM12","DZM12"))
##显示condition设置
condition
##[1] WT     WT     WT     WRKY16 WRKY16 WRKY16
##Levels: WT WRKY16
#设置colData值
#colData <- data.frame(row.names = colnames(exprSetTPM2), condition)
colData <- data.frame(row.names = colnames(exprSetTPM3), condition)
##显示colData值
colData
#condition
#WT_rep1                   WT
#WT_rep2                   WT
#WT_rep3                   WT
#GhWRKY16.RNAi_rep2    WRKY16
#GhWRKY16.RNAi_rep3    WRKY16
#GhWRKY16.RNAi_rep5    WRKY16
#构建dds矩阵
dds <- DESeqDataSetFromMatrix(exprSetTPM3, colData, design = ~condition)
#对原始dds进行normalize
dds <- DESeq(dds)
#显示dds信息
dds
##使用DESeq2包中的results()函数，提取差异分析的结果
##Usage:results(object, contrast, name, .....）
##将提取的差异分析结果定义为变量"res" 
##contrast: 定义谁和谁比较
#res = results(dds, contrast=c("condition", "CK", "CK1"))
#res = results(dds, contrast=c("condition", "SK", "SK1"))
#res = results(dds, contrast=c("condition", "S", "S1"))
#res = results(dds, contrast=c("condition", "H15", "DH15"))
res = results(dds, contrast=c("condition", "ZM12", "DZM12"))
##对结果res利用order()函数按pvalue值进行排序
##创建矩阵时，X[i,]指矩阵X中的第i行，X[,j]指矩阵X中的第j列
#order()函数先对数值排序，然后返回排序后各数值的索引，常用用法：V[order(V)]或者df[order(df$variable),]
res = res[order(res$pvalue),]
#显示res结果首信息
head(res)
#log2 fold change (MLE): condition WT vs WRKY16 
#Wald test p-value: condition WT vs WRKY16 
#DataFrame with 6 rows and 6 columns
#baseMean log2FoldChange     lfcSE      stat       pvalue         padj
#<numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
 # Gohir.A01G010800  1221.149       -6.30991  0.224838  -28.0642 2.67779e-173 1.16109e-168
#Gohir.A01G009380   360.635        4.14413  0.228304   18.1518  1.24219e-73  2.69306e-69

#对res矩阵进行总结，利用summary命令统计显示一共多少个genes上调和下调
summary(res)
#out of 59243 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 370, 0.62%
#LFC < 0 (down)     : 598, 1%
#outliers [1]       : 33, 0.056%
#low counts [2]     : 15850, 27%
#(mean count < 6)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#将分析的所有结果进行输出保存
#write.csv(res, file="All_20dparesults.csv")
#显示显著差异的数目  
table(res$padj<0.05)
#FALSE  TRUE 
#342684   676 
#使用subset()函数过滤需要的结果至新的变量diff_gene_Group2中
#根据不同参数筛选差异表达基因（需要自己调整）
diff_gene_Group2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
#diff_gene_Group2 <- subset(res, padj < 0.05 )
#也可以将差异倍数分开来写：
#> diff_gene_Group2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
#使用dim函数查看该结果的维度、规模
dim(diff_gene_Group2)
#[1] 490   6
#显示结果的首信息
summary(diff_gene_Group2)
head(diff_gene_Group2)
#log2 fold change (MLE): condition WT vs WRKY16 
#Wald test p-value: condition WT vs WRKY16 
#DataFrame with 6 rows and 6 columns
#baseMean log2FoldChange     lfcSE      stat       pvalue         padj
#<numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#  Gohir.A01G010800  1221.149       -6.30991  0.224838  -28.0642 2.67779e-173 1.16109e-168
#Gohir.A01G009380   360.635        4.14413  0.228304   18.1518  1.24219e-73  2.69306e-69
#将结果进行输出保存
write.csv(diff_gene_Group2, file = "SK6_SK6.csv")

