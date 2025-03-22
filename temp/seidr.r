##RStudio  获取去批次后TPM和基因ID
rm(list=ls())
library(data.table)
expr1=fread("C:/Users/hr345/Desktop/all/alltpm1.txt",header=T,sep="\t")
sample=fread("C:/Users/hr345/Desktop/all/WGCNA/抗感/抗感sample1.txt",header=T,sep="\t")
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
write.table(ad.ComBat,"WGCNAtpm.txt",sep = "\t",row.names = TRUE, col.names = TRUE, quote = TRUE)

##Linux
#基因ID文件名为geneid.txt，TPM文件名为WGCNAtpm.txt(格式为：行名是处理即SSR，列名是基因ID)
#删除文本的第一行内容
	sed -i "1d" WGCNAtpm.txt
#删除第一列
	cut -c 2- WGCNAtpm.txt > finalTPM.txt
#查看文件格式，如果显示“ASCII text，with CRLF line terminators”则需要转换格式，否则会报错
	file <文件名>
#改成UNIX格式
	dos2unix <文件名>
#最终格式详见: seidr相应文件格式.txt Figure 1 & Figure 2
	

##构建整体网络
#分别用13种算法进行单独计算
	vi seidr.sh   #(每个单独建立一个脚本跑比较快，最终输出的是.tsv文件)
	#一定要加绝对路径(which correlation)
	#文件名前也要加绝对路径，否则他会自动定位到/public/
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/correlation -m pearson -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/correlation -m spearman -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/pcor -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/mi -m RAW -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/mi_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/mi -m CLR -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt  -M /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/mi_scores.tsv -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/clr_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/mi -m ARACNE -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -M /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/mi_scores.tsv -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/aracne_scores.tsv
	#前面几个比较快，后面的比较慢
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/narromi -m interior-point -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/narromi_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/plsnet -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/plsnet_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/Rspace/bin/llr-ensemble -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/llr_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/svm-ensemble -k POLY -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/svm_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/genie3 -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/genie3_scores.tsv --scale
	#下面两种算法速度太慢，感觉可以舍弃
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/tigress -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/tigress_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/el-ensemble -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/elnet_scores.tsv --scale
	#.tsv文件格式详见seidr相应文件格式.txt Figure 3
	
#对这13种算法计算出的网络进行排名
	vi import.sh#(每个单独建立一个脚本跑比较快，最终输出的是.sf文件)
	#!/bin/bash
	#PBS -N import            
	#PBS -l nodes=1:ppn=2         
	#PBS -l walltime=2400:00:00 
	#PBS -l mem=20G             
	#PBS -q batch                 
	#PBS -V                      
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -A -r -u -n PEARSON -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/person_scores.sf -F lm -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/pearson_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -A -r -u -n SPEARMAN -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/spearman_scores.sf -F lm -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/spearman_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -A -r -u -n PCOR -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/pcor_scores.sf -F lm -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/pcor_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -u -n MI -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/mi_scores.sf -F lm -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/mi_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -u -z -n CLR -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/clr_scores.sf -F lm -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/clr_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -u -z -n ARACNE -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/aracne_scores.sf -F lm -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/aracne_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n NARROMI -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/narromi_scores.sf -F m -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/narromi_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n PLSNET -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/plsnet_scores.sf -F m -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/plsnet_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n LLR -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/llr_scores.sf -F m -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/llr_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n SVM -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/svm_scores.sf -F m -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/svm_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n GENIE3 -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/genie3_scores.sf -F m -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/genie3_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	#.sf文件格式详见seidr相应文件格式.txt Figure 4
	#提交脚本时用“qsub -q smp import.sh”，如果用low节点会报错

#聚合13种算法（得到aggregated.sf）
	vi aggregate.sh
	#!/bin/bash
	#PBS -N aggregate            
	#PBS -l nodes=1:ppn=2         
	#PBS -l walltime=2400:00:00 
	#PBS -l mem=20G             
	#PBS -q batch                 
	#PBS -V 
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr aggregate -m irp -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/aracne_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/clr_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/genie3_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/llr_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/mi_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/narromi_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/pcor_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/person_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/plsnet_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/spearman_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/svm_scores.sf
	#具体文件格式详见seidr相应文件格式.txt Figure 5

#修建网络
	vi backbone.sh
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr backbone -F 1.28 /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated_gdh2.sf



##构建某一目标基因网络
#文件名之前最好加上绝对路径，否则可能会定位到public路径导致报错

#选择目标基因
	echo "Gohir.D03G104800" > GDH2.txt
 
#Inferring sub-networks(用不同推断方法计算每个基因与目标基因的关联度，构建不同的子网络)，生成.tsv文件。最好每个算法单独同时跑。
	vi Seidr.sh  
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/correlation -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -m pearson -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt --scale -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/pearson_gdh2_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/correlation -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -m spearman -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/spearman_gdh2_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/pcor -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt --scale -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/pcor_gdh2_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/mi -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -m RAW -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -M /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_full_scores.tsv -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_gdh2_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/mi -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -m CLR -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt  -M /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_full_scores.tsv -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/clr_gdh2_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/mi -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -m ARACNE -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -M /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_full_scores.tsv -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aracne_fzf1_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/narromi -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -m interior-point -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/narromi_gdh2_scores.tsv
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/plsnet -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/plsnet_gdh2_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/llr-ensemble -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/llr_gdh2_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/svm-ensemble -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -k POLY -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/svm_gdh2_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/genie3 -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/genie3_gdh2_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/tigress -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/tigress_gdh2_scores.tsv --scale
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/el-ensemble -t /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/GDH2.txt -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/finalTPM.txt -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/elnet_gdh2_scores.tsv --scale
 
#importing(速度比前一步要快很多)，生成.sf文件
	vi import.sh
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -A -r -u -n PEARSON -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/person_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/pearson_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -A -r -u -n SPEARMAN -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/spearman_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/spearman_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -A -r -u -n PCOR -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/pcor_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/pcor_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -u -n MI -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -u -z -n CLR -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/clr_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/clr_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -u -z -n ARACNE -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aracne_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aracne_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -u -n MI -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -u -z -n CLR -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/clr_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/clr_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n NARROMI  -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/narromi_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/narromi_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n PLSNET -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/plsnet_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/plsnet_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n LLR -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/llr_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/llr_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n SVM -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/svm_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/svm_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n GENIE3 -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/genie3_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/genie3_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n TIGRESS -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/tigress_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/tigress_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr import -r -z -n ELNET -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/elnet_gdh2_scores.sf -F el -i /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/elnet_gdh2_scores.tsv -g /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/geneid.txt
 
#Aggregating(将第一步构建的子网络聚合成一个网络)，生成“aggregated_gdh2.sf”文件(不写入.sh会killed)
	vi aggregate.sh
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr aggregate -m irp -o /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated_gdh2.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aracne_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/clr_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/elnet_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/genie3_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/llr_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/mi_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/narromi_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/pcor_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/person_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/plsnet_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/spearman_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/svm_gdh2_scores.sf /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/tigress_gdh2_scores.sf
 
#可以查看与目标基因相关的前几个基因（最后一行出现的基因与目标基因的关联度最高，最后一列是13种方法聚合之后的综合权重和排名，所以后面作图之类的用最后一列就可以）
	seidr top -n 3 aggregated_gdh2.sf
	#Gohir.D03G104800	Gohir.A08G177600	Directed	0.699017;2	6.71211;2	0.514;2	0.40955;7      0.472;50	0.699017;2	nan;nan	0.000301451;694	0.923401;5	0.00398083;2129	0.941733;2	0.476;54.5	0.3816;0.948903;3
	#Gohir.D03G104800	Gohir.A08G213100	Directed	nan;nan	6.54461;8	0.468;13	0.240768;50    0.523;2	0.642258;13	0.472689;2	0.000446945;73	0.911057;8	0.0143462;18	0.921586;13	0.514;1	0.15275;11	0.972216;2
	#Gohir.D03G104800	Gohir.A02G159200	Directed	0.699322;1	7.17664;1	0.515;1	0.0948764;297  0.492;17	0.699322;1	nan;nan	0.000381268;197	0.94602;1	0.00989551;90	0.945725;1	0.488;30.5	0.93065;1	1;1

#Backboning(修剪网络)(不写入.sh会killed)(生成“aggregated_gdh2.bb.sf”)
	vi backbone.sh
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr backbone -F 1.28 /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated_gdh2.sf

#可以查看网络中的边缘
	seidr view --column-headers aggregated_gdh2.bb.sf | head -n 3 | column -t
	#Source            Target          Type      ARACNE_score;ARACNE_rank  CLR_score;CLR_rank  ELNET_score;ELNET_rank  GENIE3_score;GENIE3_rank  LLR_score;LLR_rank  MI_score;MI_rank  NARROMI_score;NARROMI_rank  PCOR_score;PCOR_rank  PEARSON_score;PEARSON_rank  PLSNET_score;PLSNET_rank  SPEARMAN_score;SPEARMAN_rank  SVM_score;SVM_rank  TIGRESS_score;TIGRESS_rank  irp_score;irp_rank    NC_Score;NC_SDev
#Gohir.D03G104800  Gohir.1Z014800    Directed          nan;nan                   nan;nan             nan;nan                 nan;nan              0.002;29440.5       0.0560465;51440    0.0211127;7982              -8.12e-05;22308       -0.203338;41856             0.000880798;37980         -0.226561;41101                   nan;nan                 nan;nan             0.132468;47642      0.333333;5.45132e-06
#Gohir.D03G104800  Gohir.1Z015600    Directed          nan;nan                   nan;nan             nan;nan                 nan;nan               nan;nan             0.057954;50773    0.0158132;10004             -4.23142e-05;36566    -0.176441;43878             0.000764945;41375         -0.155673;46012                   nan;nan                 nan;nan             0.0665913;55419     0.333333;3.86506e-06

 #查看特定的节点和边缘
	vi index.sh  #(创建索引，必须要进行的一步)(生成“aggregated_gdh2.bb.sf.sfi”)
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr index /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated_gdh2.bb.sf
	vi view.sh  #(还有些问题)
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr view -n Gohir.1Z014800 /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated_gdh2.bb.sf
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr view -n Gohir.1Z014800:Gohir.1Z015600 /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated_gdh2.bb.sf

#图形和中心度统计
	seidr reheader aggregated_gdh2.bb.sf  ##将断开的节点删除
	vi graph.sh    ##计算图形总体表现
	/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr graphstats /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated_gdh2.bb.sf
	##在“graph.sh.o539400”文件中可以看到具体信息
	###Number of Nodes:	56762
	###Number of Edges:	56761
	###Number of Connected Components:	1
	###Global clustering coefficient:	0
	###Scale free fit:	0.348233
	###Average degree:	1.99996
	###Average weighted degree:	0.426842
	###Network diameter:	1.97222
	###Average path length:	0.426835
	
#计算节点中心性
	vi centrality.sh
	#(存在问题)/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/seidr stats --exact /vol3/agis/huguanjing_group/xiongxianpeng/HR_transcriptome/seidr/GDH2/aggregated_gdh2.bb.sf
	seidr view --centrality aggregated_gdh2.bb.sf | sort -k2g | tail -n 5 | column -t   ##查看数据
