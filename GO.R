rm(list=ls())
setwd("C:/Users/hr345/Desktop/all/GO")
library(clusterProfiler)
library(stringr)
library(dplyr)
library(tidyverse)
library(stringr)
library(KEGGREST)
library(AnnotationForge)
#egg_f <- "Ghirsutum_527_v2.1.protein_xxp.annotations"
#egg1 <- read.csv(egg_f, sep = "\t")
#dim(egg1)
#[1] 101038     22
#gene56763<-read.table("56763.txt", sep = "\t",header = T)
#head(gene56763)
#query_name
#1 Gohir.1Z001700.1
#egg56763=merge(gene56763,egg1,by="query_name")
#egg<-egg56763
egg_f <- "egg.annotations"
egg1 <- read.csv(egg_f,sep = "\t")
gene56763<-read.table("56763_2.txt", sep = "\t",header = T)
egg56763=merge(gene56763,egg1,by="query_name")
egg1<-egg56763
dim(egg1)
#
gene_ids <- egg1$query_name
eggnog_lines_with_go <- egg1$GOs!= ""
eggnog_lines_with_go
eggnog_annoations_go <- str_split(egg1[eggnog_lines_with_go,]$GOs, ",")
gene_to_go <- data.frame(gene = rep(gene_ids[eggnog_lines_with_go],
                                    times = sapply(eggnog_annoations_go, length)),
                         term = unlist(eggnog_annoations_go))
head(gene_to_go)
# gene       term
#1 Gohir.1Z015700.1 GO:0005575
#2 Gohir.1Z015700.1 GO:0005622
dim(gene_to_go)
#[1] 2035409       2
#save(gene_to_go,file = "GO.RData")
#load(file = "GO.RData")



gene_list<-read.table("C:/Users/hr345/Desktop/all/GO/turquoise/turquoise.txt", sep = "\t")
gene_list<-read.table("C:/Users/hr345/Desktop/all/GO/blue/blue.txt", sep = "\t")
gene_list<-read.table("C:/Users/hr345/Desktop/all/GO/black/black.txt", sep = "\t")
gene_list<-read.table("C:/Users/hr345/Desktop/all/GO/purple/purple.txt", sep = "\t")
gene_list<-read.table("C:/Users/hr345/Desktop/all/GO/brown/brown.txt", sep = "\t")
gene_list=gene_list$V1
term2gene<-gene_to_go[,c(2,1)]
df<-enricher(gene=gene_list,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             TERM2GENE = term2gene)
head(df)
barplot(df)
dotplot(df)
df<-as.data.frame(df)
head(df)
dim(df)
df1<-go2term(df$ID)
dim(df1)
head(df1)
class(df1)
names(df1) <- c("ID","term")
dff=merge(df1,df,by="ID")

dff$term<-df1$term
df2<-go2ont(dff$ID)
dim(df2)
head(df2)
dff$Ont<-df2$Ontology
head(dff)
write.table(dff,"C:/Users/hr345/Desktop/all/GO/turquoise/turquoise2.GO.txt",sep = "\t")
write.table(dff,"C:/Users/hr345/Desktop/allGO/blue/blue2.GO.txt",sep = "\t")
write.table(dff,"C:/Users/hr345/Desktop/all/GO/black/black2.GO.txt",sep = "\t")
write.table(dff,"C:/Users/hr345/Desktop/all/GO/purple/purple2.GO.txt",sep = "\t")
write.table(dff,"C:/Users/hr345/Desktop/all/GO/brown/brown2.GO.txt",sep = "\t")

##pheatmap
aa=fread("C:/Users/hr345/Desktop/all/GO/GO results.txt",header=T,sep="\t")
aa <- as.data.frame(aa)
class(aa)
rownames(aa)=aa$term
aa=aa[,-1]
library(ggplot2)
library(pheatmap)
pheatmap(aa,scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,fontsize_col = 10,fontsize_row = 8,show_colnames = T,cellwidth=23,cellheight = 7,angle_col = 45,color = colorRampPalette(c("white","red"))(100))

library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_scout()
df3<-dff%>%
  select(c("term","Ont","pvalue"))
head(df3)
library(ggplot2)
ggplot(df3[1:20,],aes(x=term,y=-log10(pvalue)))+
  geom_col(aes(fill=Ont))+
  coord_flip()+labs(x="")+
  theme_bw()
