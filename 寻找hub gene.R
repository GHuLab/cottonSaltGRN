###寻找hub gene
rm(list=ls())
library(RColorBrewer)
library(flashClust)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
load("C:/Users/hr345/Desktop/all/WGCNA/抗感/Rtemp56763.rdata")
load("C:/Users/hr345/Desktop/all/WGCNA/抗感/wgcna.56763.cgn.Rdata")
load("C:/Users/hr345/Desktop/all/WGCNA/抗感/wgcna.TPM.56763.prep.Rdata")

ls()
i="cgnP22b"
net = get(i)
datExpr=as.data.frame(t(expr))
dim(datExpr)
# [[1]   123 56763
moduleLabels = net$colors

moduleColors = labels2colors(net$colors)
table(moduleColors)

##moduleColors
#       black         blue        brown         cyan        green 
#         347         4540         4264          138          485 
# greenyellow         grey    lightcyan      magenta midnightblue 
#         237        27128          120          258          130 
#        pink       purple          red       salmon          tan 
#         300          253          387          148          173 
#   turquoise       yellow 
#       16588         1267 

MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(MEs)
dim(MET)
TOM = TOMsimilarityFromExpr(datExpr, power = 22)
module = "brown"
probes = names(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""), nodeFile = paste("CytoscapeInput-nodes-", paste(module,  collapse="-"), ".txt", sep=""),weighted = TRUE, threshold = 0.1,nodeNames = modProbes, nodeAttr = moduleColors[inModule])
##KME
modNames = substring(names(MEs), 3)
datKME=signedKME(datExpr, MEs, outputColumnName="kME_MM.")
module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module
blue_module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes])
names(blue_module)="genename"
blue_KME<-as.data.frame(datKME[moduleGenes,column])
names(blue_KME)="KME"
rownames(blue_KME)=blue_module$genename
FilterGenes = abs(blue_KME$KME) > 0.8
table(FilterGenes)
blue_hub<-subset(blue_KME, abs(blue_KME$KME)>0.8)
write.csv(blue_hub, "hubgene_KME_blue.csv")

#不用看
##脂肪酸通路基因表达量热图
rm(list=ls())
MEs2=fread("C:/Users/hr345/Desktop/脂肪酸表达量1.txt",header=T,sep="\t")
MEs2=MEs2[,-1]
c=data.frame(cbind(apply(MEs2[,1:13],1,mean),
                   apply(MEs2[,14:16],1,mean),
                   apply(MEs2[,17:19],1,mean),
                   apply(MEs2[,20:22],1,mean),
                   apply(MEs2[,23:25],1,mean),
                   apply(MEs2[,26:28],1,mean),
                   apply(MEs2[,29:31],1,mean),
                   apply(MEs2[,32:34],1,mean),
                   apply(MEs2[,35:37],1,mean),
                   apply(MEs2[,38:40],1,mean),
                   apply(MEs2[,41:43],1,mean),
                   apply(MEs2[,44:46],1,mean),
                   apply(MEs2[,47:49],1,mean),
                   apply(MEs2[,50:52],1,mean),
                   apply(MEs2[,53:55],1,mean),
                   apply(MEs2[,56:58],1,mean),
                   apply(MEs2[,59:61],1,mean),
                   apply(MEs2[,62:64],1,mean),
                   apply(MEs2[,65:67],1,mean),
                   apply(MEs2[,68:70],1,mean),
                   apply(MEs2[,71:73],1,mean),
                   apply(MEs2[,74:76],1,mean),
                   apply(MEs2[,77:79],1,mean),
                   apply(MEs2[,80:82],1,mean),
                   apply(MEs2[,83:85],1,mean),
                   apply(MEs2[,86:88],1,mean),
                   apply(MEs2[,89:91],1,mean),
                   apply(MEs2[,92:94],1,mean),
                   apply(MEs2[,95:96],1,mean),
                   apply(MEs2[,97:98],1,mean),
                   apply(MEs2[,99:101],1,mean),
                   apply(MEs2[,102:104],1,mean),
                   apply(MEs2[,105:108],1,mean),
                   apply(MEs2[,109:111],1,mean),
                   apply(MEs2[,112:114],1,mean),
                   apply(MEs2[,115:117],1,mean),
                   apply(MEs2[,118:120],1,mean),
                   apply(MEs2[,121:123],1,mean)))

class(c)
#[1] "data.frame"
#修改行名
names(c) <- c("a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","p1","q1","r1","s1")
c$module=c("Gohir.D03G060500(SAD6)","Gohir.A02G107100(SAD6)","Gohir.D13G172400","Gohir.A13G161600","Gohir.A11G244900","Gohir.D07G003300(FATB)","Gohir.A07G003000(FATB)","Gohir.A06G060100(FATB)","Gohir.D06G058300(FATB)", "Gohir.D08G212700(FATA)","Gohir.A08G195400(FATA)","Gohir.D06G045200(LACS9)","Gohir.A06G044700(LACS9)","Gohir.A05G234900(LACS9)","Gohir.A07G087200(LACS9)", "Gohir.D09G089900(FAD3)","Gohir.D11G331500(FAD2)","Gohir.A11G314400(FAD2)","Gohir.A01G137630(FAD2)","Gohir.D01G125900(FAD2)","Gohir.A06G201800(LPCAT)","Gohir.D06G220000(LPCAT)","Gohir.D02G105800(LPCAT)","Gohir.A03G086500(LPCAT)","Gohir.D04G170000(GPAT9)","Gohir.A04G125200(GPAT9)","Gohir.A03G209200(GPAT9)","Gohir.A01G043500(LPP2)","Gohir.D01G034600(LPP2)","Gohir.A09G031100(GPDH)","Gohir.D09G031300(GPDH)","Gohir.A13G189400(KAS2)","Gohir.D13G196305(KAS2)","Gohir.A08G246500(KAS2)","Gohir.D08G267700(KAS2)","Gohir.D11G059800(DGAT2)","Gohir.A11G055700(DGAT2)","Gohir.D11G1092000(DGAT3)","Gohir.A11G104400(DGAT3)")
c <- c[, c("module","a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","p1","q1","r1","s1")]

anno=read.table("C:/Users/hr345/Desktop/all/WGCNA/抗感/抗感分类.txt",header=T,sep="\t")
Period=anno[,-1]
#names(Period)="time"
head(Period)
#Period  DPA
#1 Elongation 5dpa
class(Period)
dim(Period)
#anno11=data.frame(Period) 非数值需要转置


per=c
dat2 <- per[,-1]*10
rownames(Period) = rownames(t(dat2))
row.names(dat2)=per$module
#dat2=t(dat2)
library(pheatmap)


Period$time<- factor(Period$time, levels=c("0h", "15m","1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
Period$treatment<- factor(Period$treatment, levels=c("mock","150(T)","200(T)","200(S)","200(TM-1)","salt-alkali(S)","salt-alkali(T)","sk(T)"), ordered=TRUE)

ann_colors1 = list(time = c("0h"="#fffff5", "15m"= "#fab1ce","1h"="#EE7785","3h"="#ef5285","6h"="#C5E99B","12h"="#8CD790","24h"="#77AF9C","48h"="#379392","72h"="#285943"),
                   treatment = c("mock"="#fc913a", "150(T)"= "#9baec8","200(T)"="#47b8e0","200(S)"="#2b90d9","200(TM-1)"="#4F86C6","salt-alkali(S)"="#dedcee","salt-alkali(T)"="#6a60a9","sk(T)"="#5c196b"))
pheatmap(dat2,annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 7,show_colnames = F,
         gaps_col = c(5,9,16,23,28,31,34,38))


#卡方检验
R=matrix(c(120,56763,9,4264),nrow = 2,ncol = 2)
chisq.test(R)
