##PCA
rm(list=ls())
library(data.table)
expr1=fread("F:/new/expr2.txt",header=T,sep="\t")
sample =fread("C:/Users/hr345/Desktop/all/new/sample1.txt",header=T,sep="\t")
colnames(expr1)[1] <- 'sample'
aa=merge(sample,expr1,by="sample",sort=F)
class(aa)
#[1] "data.table" "data.frame"
aa <- as.data.frame(aa)
rownames(aa)=aa$sample
expr=t(aa[,-1:-3])
dim(expr)
#[1] 75376    18
anno=sample
dim(anno)
#[1] 18  3
expr1=t(expr)
expr=t(expr1)
dim(expr)
#[1] 75376    18
expr1=t(expr)
b=apply(expr1,2,as.numeric)
df_pca <- prcomp(log2(b+1))
anno$treatment<- factor(anno$treatment, levels=c("A", "S","CK"), ordered=TRUE)
anno$TRV<- factor(anno$TRV, levels=c("GDH2","mock"), ordered=TRUE)
df_pcs <-data.frame(df_pca$x,treatment = anno$treatment,TRV=anno$TRV) 
head(df_pcs,3)
percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)
percentage<-paste(colnames(df_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(df_pcs,aes(x=PC1,y=PC2,color=treatment,shape=TRV))+
     geom_point(size=3)+ 
     xlab(percentage[1]) +
     ylab(percentage[2])+
     theme(plot.title = element_text(size = 1, face = "bold") ,legend.title=element_text(size=15), 
           legend.text=element_text(size=15)) +theme(axis.title.x = element_text(size = 15)) +
     theme(axis.title.y = element_text(size = 15))+ 
     scale_color_manual(values=c("#D1B6E1","#4ea1d3","#C5E99B"))+ scale_shape_manual(values=c(15,16))

#乙烯通路热图
c <- read.table(file = "F:/new/新建文本文档.txt", sep = "\t", header = T, row.names= 1 )
c=data.frame(cbind(apply(c[,1:3],1,mean),
                    apply(c[,4:6],1,mean),
                    apply(c[,7:9],1,mean),
                    apply(c[,10:12],1,mean),
                    apply(c[,13:15],1,mean),
                    apply(c[,16:18],1,mean)))
c <- read.table(file = "D:/乙烯average.txt", sep = "\t", header = T, row.names= 1 )
pheatmap(c,scale="row",cluster_cols=F,cluster_rows =T,display_numbers = F,fontsize_col = 10,fontsize_row = 8,show_colnames = T,cellwidth=20,cellheight = 10,angle_col = 45,color = colorRampPalette(c("blue","white","red"))(10))

#火山图
rm(list=ls())
library(xlsx)
data <- read.xlsx("F:/new/CK差异表达基因.xlsx",1,encoding="UTF-8")
rownames(data)=data[,1]
data=data[,-1]
cut_off_padj =0.05
cut_off_log2FoldChange =1
data$Sig = ifelse(data$padj < cut_off_padj &
                  abs(data$log2FoldChange) >= cut_off_log2FoldChange,
                  ifelse(data$log2FoldChange > cut_off_log2FoldChange ,'Up','Down'),'no')
data = data.frame(data)
table(data$Sig)
#Down   Up 
#595   2151
library(ggplot2)
pl <- ggplot(data, aes(x =log2FoldChange, y=padj_2, colour=Sig)) +
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#45d9fd", "#ff5f2e")) + xlim(c(-10, 10)) +  #调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FoldChange,cut_off_log2FoldChange),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_padj), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="log2FoldChange", y="-log10FDR") +  #x、y轴标签
  ggtitle("TRV:00 : TRV:GhGDH2") + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right", 
        legend.title = element_blank())
pl

##GDH沉默效率柱状图
library(ggprism)
library(ggplot2)
library(data.table)
library(ggpubr)
expr1=fread("F:/new/new/GDH沉默效率.txt",header=T,sep="\t")
ggbarplot(expr1,'group','mrna',fill = "black",add = "mean_sd",xlab = F,ylab = 'Relative mRNA expression',legend='none',ggtheme = theme_prism())+
stat_compare_means(aes(label = ..p.signif..),  ## 改成星星
comparisons = list(c('TRV:00','TRV:GhGDH2')),  ## 添加一下列表
method = 't.test')

#MDA
rm(list=ls())
library(ggprism)
library(ggplot2)
library(data.table)
library(ggpubr)
expr1=fread("F:/new/new/生理指标/MDA.txt",header=T,sep="\t")
ggbarplot(expr1,x="treatment",y="MDA",color="group",fill="group",add = "mean_sd",xlab = F,ylab = "MDA cotent",ggtheme = theme_prism(), position = position_dodge(0.8),palette = c("grey","black"))+
stat_compare_means(aes(group=group),method = "t.test")
#叶绿素a
rm(list=ls())
expr1=fread("F:/new/new/生理指标/叶绿素a.txt",header=T,sep="\t")
ggbarplot(expr1,x="treatment",y="chlorophyll A",color="group",fill="group",add = "mean_sd",xlab = F,ylab = "leaf chlorophyll-A content",ggtheme = theme_prism(), position = position_dodge(0.8),palette = c("grey","black"))+
  stat_compare_means(aes(group=group),method = "t.test",label = "p.signif")
#叶绿素B
rm(list=ls())
expr1=fread("F:/new/new/生理指标/叶绿素b.txt",header=T,sep="\t")
ggbarplot(expr1,x="treatment",y="chlorophyll b",color="group",fill="group",add = "mean_sd",xlab = F,ylab = "leaf chlorophyll-b content",ggtheme = theme_prism(), position = position_dodge(0.8),palette = c("grey","black"))+
  stat_compare_means(aes(group=group),method = "t.test",label = "p.signif")
#总叶绿素
rm(list=ls())
expr1=fread("F:/new/new/生理指标/总叶绿素.txt",header=T,sep="\t")
ggbarplot(expr1,x="treatment",y="chlorophyll",color="group",fill="group",add = "mean_sd",xlab = F,ylab = "leaf chlorophyll content",ggtheme = theme_prism(), position = position_dodge(0.8),palette = c("grey","black"))+
  stat_compare_means(aes(group=group),method = "t.test",label = "p.signif")
#POD
rm(list=ls())
expr1=fread("F:/new/new/生理指标/POD.txt",header=T,sep="\t")


#带有误差线的折线图
rm(list=ls())
library(ggprism)
library(ggplot2)
library(data.table)
library(ggpubr)
expr1=fread("F:/all1/GDH2表达量/200.txt",header=T,sep="\t")
ggline(expr1,x = "time",y = "expr",add = "mean_se",color="species",xlab = "time",ylab = "GhGDH2 expression")
