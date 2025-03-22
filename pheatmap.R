##salt-alkali
rm(list=ls())
c <- read.table(file = "C:/Users/hr345/Desktop/all/GDH expression/salt-alkali/GDH-average.txt", sep = "\t", header = T, row.names= 1 )
class(c)
names(c) <- c("a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","p1","q1","r1")
c <- c[, c("a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","p1","q1","r1")]
anno=read.table("C:/Users/hr345/Desktop/all/GDH expression/salt-alkali/sample.txt",header=T,sep="\t")
Period=anno[,-1]
head(Period)
class(Period)
dim(Period)
per=c
rownames(Period) = rownames(t(per))
library(pheatmap)
Period$time<- factor(Period$time, levels=c("15m","1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
ann_colors1 = list(time = c( "15m"= "#fab1ce","1h"="#EE7785","3h"="#ef5285","6h"="#C5E99B","12h"="#8CD790","24h"="#77AF9C","48h"="#379392","72h"="#285943"),treatment = c("mock"="#fc913a", "150(T)"= "#9baec8","200(T)"="#47b8e0","200(S)"="#2b90d9","200(TM-1)"="#4F86C6","salt-alkali(S)"="#dedcee","salt-alkali(T)"="#6a60a9","sk(T)"="#5c196b"))

pheatmap(per,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,show_colnames = F,cellwidth=15,cellheight = 30,
         gaps_col = c(4,8,15,22,27,30,33,37),
         color = colorRampPalette(c("blue","white","red"))(100))

pheatmap(per,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="column",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,show_colnames = F,cellwidth=15,cellheight = 30,
         gaps_col = c(4,8,15,22,27,30,33,37),
         color = colorRampPalette(c("blue","white","red"))(100))

pheatmap(per,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="row",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,show_colnames = F,cellwidth=15,cellheight = 30,
         gaps_col = c(4,8,15,22,27,30,33,37),
         color = colorRampPalette(c("blue","white","red"))(100))

##different treatment
rm(list=ls())
c <- read.table(file = "C:/Users/hr345/Desktop/all/GDH expression/different treatment/all.txt", sep = "\t", header = T, row.names= 1 )
names(c) <- c("a", "b", "c","d","e","f","g","h","i","j","k","a1", "b1", "c1","d1","e1","f1","g1","h1","i1")
c <- c[, c("a", "b", "c","d","e","f","g","h","i","j","k","a1", "b1", "c1","d1","e1","f1","g1","h1","i1")]
anno=read.table("C:/Users/hr345/Desktop/all/GDH expression/different treatment/sample1.txt",header=T,sep="\t")
Period=anno[,-1]
head(Period)
dim(Period)
per=c
rownames(Period) = rownames(t(per))
library(pheatmap)
Period$time<- factor(Period$time, levels=c( "1h","3h","6h","12h","24h"), ordered=TRUE)
ann_colors1 = list(time = c("1h"="#EE7785","3h"="#ef5285","6h"="#C5E99B","12h"="#8CD790","24h"="#77AF9C"),treatment = c("mock"="#fc913a", "PEG"= "#9baec8","37C"="#47b8e0","4C"="#6a60a9"))

pheatmap(per,annotation_col = Period,
          annotation_colors = ann_colors1,
          scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
          fontsize = 15,show_colnames = F,cellwidth=20,cellheight = 30,
          gaps_col = c(5,10,15,20),color = colorRampPalette(c("blue","white","red"))(100))

pheatmap(per,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="row",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,show_colnames = F,cellwidth=20,cellheight = 30,
         gaps_col = c(5,10,15,20),
         color = colorRampPalette(c("blue","white","red"))(100))

pheatmap(per,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="column",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,show_colnames = F,cellwidth=15,cellheight = 30,
         gaps_col = c(5,10,15,20),color = colorRampPalette(c("blue","white","red"))(100))




##different tissue
rm(list=ls())
c <- read.table(file = "C:/Users/hr345/Desktop/all/GDH expression/different tissue/new.txt", sep = "\t", header = T, row.names= 1 )
names(c) <- c("a", "b", "c","d","e","f","g","h","i","a1", "b1", "c1","d1","e1","f1","g1","h1","i1")
c <- c[, c("a", "b", "c","d","e","f","g","h","i","a1", "b1", "c1","d1","e1","f1","g1","h1","i1")]
anno=read.table("C:/Users/hr345/Desktop/all/GDH expression/different tissue/newsample.txt",header=T,sep="\t")
Period=anno[,-1]
head(Period)
class(Period)
per=c
rownames(Period) = rownames(t(per))
Period$tissue<- factor(Period$tissue, levels=c("root", "stem","leaf","bract","sepal","torus"), ordered=TRUE)
ann_colors1 = list(treatment = c("root"="#D499B9", "stem"= "#9055A2","leaf"="#791E94","bract"="#8CD790","sepal"="#77AF9C","torus"="#379392"))
pheatmap(per,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="row",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,show_colnames = F,cellwidth=15,cellheight = 30,
         gaps_col = c(3,6,9,12,15,18),color = colorRampPalette(c("midnightblue","white","red"))(10))
