
########WGCNA

rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/WGCNA/sample1.txt",header=T,sep="\t")
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

srr12=anno
coldata1<-data.frame(srr=srr12$SRR,sample = srr12$SRR)
coldata=coldata1
head(coldata)
#srr      sample
#1 SRR10518033 SRR10518033
#2 SRR10518014 SRR10518014
dim(coldata)
#123 2
library(RColorBrewer)
library(WGCNA)
library(flashClust)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

multiExpr = vector(mode = "list", length = 1)
multiExpr[[1]] = list(data = t(expr))


checkSets(multiExpr)
# $nSets [1] 1
# $nGenes [1]56763
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
# $nGenes [1] 56763
# $nSamples [1] 123
# $structureOK [1] TRUE
shortLabels="All"
save(expr, coldata, multiExpr,  file="/vol3/agis/huguanjing_group/xiongxianpeng/WGCNA/Rtemp56763.rdata")

## Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
# types<-c("unsigned", "signed", "signed hybrid")
powers = c(c(1:10), seq(from = 12, to=40, by=2))
for(type in c("unsigned", "signed", "signed hybrid")){
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
for(type in c("signed", "signed hybrid"))
{
  for(corU in c("cor","bicor")){
    p=paste0("powerTables.",gsub(" ","",type),".",corU)
    print(p)
    print(lapply(get(p), function(x){ x<-x$data; x$Power[x$SFT.R.sq>0.8 & x$slope<0][1]} ))
    powerTables[[paste0(gsub(" ","",type),".",corU)]] <-get(p)[[1]]
  }
}

powers=c(10,22,24)
save(expr, multiExpr, powerTables, powers, file = "/vol3/agis/huguanjing_group/xiongxianpeng/WGCNA/wgcna.TPM.56763.prep.Rdata")


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
    save(list=grep("cgnP",ls(),value=TRUE), file = "/vol3/agis/huguanjing_group/xiongxianpeng/WGCNA/wgcna.cgn56763.Rdata")
  }
}




# 34672 genes 
library(RColorBrewer)
library(flashClust)
library(WGCNA)
load("C:/Users/hr345/Desktop/all/WGCNA/Rtemp56763.rdata")
load("C:/Users/hr345/Desktop/all/WGCNA/wgcna.56763.cgn.Rdata")
load("C:/Users/hr345/Desktop/all/WGCNA/wgcna.TPM.56763.prep.Rdata")

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
#123 14
pdf(file = "Eigengene adjacency heatmap.pdf", width = 12, height = 9)
plotEigengeneNetworks(MET, "C:/Users/hr345/Desktop/all/WGCNA/Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()
library(stringr)
MEs=net$MEs
MEs_col=MEs
colnames(MEs_col)=paste0('ME',labels2colors(
  as.numeric(str_replace_all(colnames(MEs),'ME',''))))
MEs_col=orderMEs(MEs_col)
plotEigengeneNetworks(MEs_col,'Eigengene adjacency heatmap',
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2),plotDendrograms = T,
                      xLabelsAngle=90)


# WGCNA hub gene heatmap
rm(list=ls())
library(WGCNA);
library(RColorBrewer)
library(flashClust);
load("C:/Users/hr345/Desktop/all/WGCNA/Rtemp56763.rdata")
load("C:/Users/hr345/Desktop/all/WGCNA/wgcna.56763.cgn.Rdata")
load("C:/Users/hr345/Desktop/all/WGCNA/wgcna.TPM.56763.prep.Rdata")

ls()
i="cgnP22b"
net = get(i)
datExpr=as.data.frame(t(expr))
dim(datExpr)
# [1]413 25441
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)

MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
#calculate MM
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
nSamples = nrow(datExpr)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

MET = orderMEs(MEs)
dim(MET)

# WGCNA hub gene heatmap(new)
c=data.frame(cbind(apply(MEs2[,2:4],1,mean),
                   MEs2[,1],
                   apply(MEs2[,5:7],1,mean),
                   apply(MEs2[,8:10],1,mean),
                   apply(MEs2[,11:13],1,mean),
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
names(c) <- c("a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","p1","q1","r1","s1","t1","u1")
c$module=c("black(347)","blue(4540)","brown(4264)","cyan(138)","green(485)","greenyellow(237)","grey(27128)","lightcyan(120)","magenta(258)", "midnightblue(130)","pink(300)","purple(253)","red(387)","salmon(148)","tan(173)", "turquoise(16588)","yellow(1267)")
c <- c[, c("module","a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","p1","q1","r1","s1","t1","u1")]
anno=read.table("C:/Users/hr345/Desktop/all/WGCNA/3.txt",header=T,sep="\t")
Period=anno[,-1]
#names(Period)="time"
head(Period)
#Period  DPA
#1 Elongation 5dpa
class(Period)
dim(Period)
per=c
dat2 <- per[,-1]*10
rownames(Period) = rownames(t(dat2))
row.names(dat2)=per$module
#dat2=t(dat2)
library(pheatmap)
Period$time<- factor(Period$time, levels=c("0h", "15m","1h","3h","6h","12h","24h","48h","72h"), ordered=TRUE)
Period$treatment<- factor(Period$treatment, levels=c("150(mock)","150(T)","200(T)","200(S)","200(TM-1)","salt-alkali(S)","salt-alkali(T)","sk(T)"), ordered=TRUE)
ann_colors1 = list(time = c("0h"="#fffff5", "15m"= "#fab1ce","1h"="#EE7785","3h"="#ef5285","6h"="#C5E99B","12h"="#8CD790","24h"="#77AF9C","48h"="#379392","72h"="#285943"),
                   +treatment = c("150(mock)"="#f6ea8c", "150(T)"= "#9baec8","200(T)"="#47b8e0","200(S)"="#2b90d9","200(TM-1)"="#4F86C6","salt-alkali(S)"="#dedcee","salt-alkali(T)"="#6a60a9","sk(T)"="#5c196b"))
pheatmap(dat2[-7,],annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 12,show_colnames = F,cellwidth=10,cellheight = 20,
         gaps_col = c(1,2,3,4,5,13,20,27,32,35,38,42),color = colorRampPalette(c("darkblue","white","red"))(10))



# WGCNA hub gene heatmap(old)
MEs1=t(MEs)
MEs2=MEs1
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
names(c) <- c("a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","p1","q1","r1","s1")

#c$module=names(MEs)
c$module=c("black(347)","blue(4540)","brown(4264)","cyan(138)","green(485)","greenyellow(237)","grey(27128)","lightcyan(120)","magenta(258)", "midnightblue(130)","pink(300)","purple(253)","red(387)","salmon(148)","tan(173)", "turquoise(16588)","yellow(1267)")



c <- c[, c("module","a", "b", "c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","a1", "b1", "c1","d1","e1","f1","g1","h1","i1","j1","k1","l1","m1","n1","o1","p1","q1","r1","s1")]

anno=read.table("C:/Users/hr345/Desktop/all/WGCNA/resistance.txt",header=T,sep="\t")
Period=anno[,-1]
#names(Period)="time"
head(Period)
#Period  DPA
#1 Elongation 5dpa
class(Period)
dim(Period)
#anno11=data.frame(Period) 


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



pheatmap(dat2[-7,],annotation_col = Period,
         annotation_colors = ann_colors1,
         scale="none",cluster_cols=F,cluster_rows =T,display_numbers = F,
         fontsize = 15,show_colnames = F,
         gaps_col = c(5,9,16,23,28,31,34,38))





moduleGenes = names(net$colors)[net$colors%in%c("4")]
write.csv(moduleGenes,"yellow.csv")
# 0:grey   1:turquoise   2:blue   3:brown   4:yellow   5:green    6:red    7:black
# 8:pink   9:magenta     10:purple    11:greenyellow   12:tan     13:salmon
# 14:cyan  15:midnightblue      16:lightcyan





#pdf
library(RColorBrewer)
display.brewer.all() 

brewer.pal(9,"YlOrBr")

display.brewer.pal(9,"YlGnBu")

rm(list=ls())
library(RColorBrewer)
library(flashClust)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
load("C:/Users/hr345/Desktop/all/WGCNA/抗感/Rtemp56763.rdata")
load("C:/Users/hr345/Desktop/all/WGCNA/抗感/wgcna.56763.cgn.Rdata")
load("C:/Users/hr345/Desktop/all/WGCNA/抗感/wgcna.TPM.56763.prep.Rdata")


##calculate kME
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