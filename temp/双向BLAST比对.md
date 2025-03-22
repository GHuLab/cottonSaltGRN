### 2023/9/26
根据文章“Identification of CDK gene family and functional analysis of CqCDK15 under drought and salt stress in quinoa”，**附表1**提供的拟南芥CDK基因进行RHB位点寻找，一方面需要在拟南芥和棉花间进行blast比对，另一方面需要将**附表1**数据与拟南芥Arport11进行比对找到一致序列。

1. 将附表1内容与拟南芥数据库比对，挑选bitscore最高序列的序列作为输出。
2. ~~下载Tair的公共数据，gff文件与fa文文件，其中gff文件只保留最长转录本~~，Tair用的基因组是TAIR10版本，11版本的注释，我觉得基因数少没必要在拆最长转录本，拟南芥转录本注释可能也算比较清楚。
3. 通过gffread，根据参考基因组获得蛋白序列，如果有CDS feature更好，得到两份蛋白序列文件
4. 分别建立blast数据库，棉花比拟南芥，拟南芥比棉花。输出限制e-10，输出前5个
5. 对两个结果进行筛选，挑选正反都有的结果作为输出


```bash
conda activate primer
# 将附表序列与所有的拟南芥蛋白序列比对
# diamond 建库
diamond makedb --in ../2023_5_29/dat/Araport11_pep_20220914 --db dat/ref/Ara11
diamond blastp --threads 8 --query dat/拟南芥CDK基因蛋白质序列.fasta --db dat/ref/Ara11.dmnd --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --max-target-seqs 1 --evalue 10e-10 --out res/CDKat.diamond
## diamond blastp输出的最高的结果，idenlen都是100
# ★另外再用原来的blastp尝试一下
blastp -num_threads 12 -query dat/拟南芥CDK基因蛋白质序列.fasta -db ../../ref/Ara11 -outfmt '6 qaccver saccver pident
length mismatch gapopen qstart qend sstart send evalue bitscore'  -num_alignments 5 -out res/CDKat.blastp -seg yes
# 结果上与diamond结果差不多，就是两者参数有些小差别。

# 获得UTX2.1的最长转录本与最长蛋白
# 拿到UTX2.1参考基因组与参考注释，去除scaffold
(maintoolkits) niu 10:29:47 ~/project/seqAlign/exp/2023_09_27/dat
$ vim chromosome.name
# 使用samtools去除
(maintoolkits) niu 10:36:07 ~/project/seqAlign/exp/2023_09_27/dat
$ while read line; do samtools faidx Ghirsutum_527_v2.0.fa $line >> UTX2.1.chromosome.fasta;done < chromosome.name
# 修改gff
(maintoolkits) niu 10:39:54 ~/project/seqAlign/exp/2023_09_27/dat
$ awk '/^[AD][0-9]+|^#/{print $0}' Ghirsutum_527_v2.1.gene_exons.gff3 > UTX2.1.gene_exon.gff
# 抽取最长转录本
(maintoolkits) niu 10:46:05 ~/project/seqAlign/exp/2023_09_27/dat
$ sh ../run/getlongestTrans_gff.sh UTX2.1.gene_exon.gff UTX2.1.longestTrans.gff
# 获得cds peptide序列
(maintoolkits) niu 10:49:51 ~/project/seqAlign/exp/2023_09_27/dat
$ gffread -g UTX2.1.chromosome.fasta -y UTX2.1.peptide.fa UTX2.1.longestTrans.gff

# 建立UTX diamond库
diamond makedb --in dat/UTX2.1.peptide.fa --db dat/ref/UTX2.1.longestTranspep
# Ara11 query UTX
diamond blastp --threads 8 --query dat/Ara11.pep.fa --db dat/ref/UTX2.1.longestTranspep.dmnd --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --max-target-seqs 5 --evalue 10e-10 --out res/AraVsUtx.diamond.blastp
# UTX query Ara11
diamond blastp --threads 8 --query dat/Ara11.pep.fa --db dat/ref/UTX2.1.longestTranspep.dmnd --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --max-target-seqs 5 --evalue 10e-10 --out res/AraVsUtx.diamond.blastp

# 合并结果，uniq -c 查看有两个的结果并提取，从之前的CDK基因转At名得到对应的Atlist，最后终提取RHB.CDK.list
awk -F "\t" '{print $2 "\t" $1}' res/UtxVsAra.diamond.blastp > res/RHB.list
cut -f1,2 res/AraVsUtx.diamond.blastp >> res/RHB.list
cat res/RHB.list | sort | uniq -c > res/RHB.uniq.list
grep -E "^ +2 " res/RHB.uniq.list > res/RHB.res.list
# 抽取CDK基因
awk -F "\t" '{print $2}' res/CDKat.diamond > tmp/CDKatID.list
# 提取RHB的CDK基因
while read line; do grep "${line}" res/RHB.res.list >> res/RHB.CDK.res.list; done < tmp/CDKatID.list
```

抽取前和抽取后gff文件feature比较
	(maintoolkits) niu 10:46:49 ~/project/seqAlign/exp/2023_09_27/dat
	$ cut -f3 UTX2.1.gene_exon.gff | sort | uniq -c
	      1 ##annot-version v2.1
	 613199 CDS
	 666351 exon
	 120317 five_prime_UTR
	  74902 gene
	      1 ##gff-version 3
	 106647 mRNA
	      1 ##species Gossypium hirsutum
	 110572 three_prime_UTR
	(maintoolkits) niu 10:47:39 ~/project/seqAlign/exp/2023_09_27/dat
	$ cut -f3 UTX2.1.longestTrans.gff | sort | uniq -c
	      1 ##annot-version v2.1
	 371639 CDS
	 394957 exon
	  71859 five_prime_UTR
	  74902 gene
	      1 ##gff-version 3
	  74902 mRNA
	      1 ##species Gossypium hirsutum
	  67307 three_prime_UTR

此外，还有通过orthofinder找到的结果，一起看一下。
```bash
# 查找符合染色体倍性的基因(忽略串联重复)
(primer) niu 09:58:39 /mnt/c/Users/HP/Desktop/张诗诗/OG结果
$ awk -F "\t" '{if(split($2,arr,", ") < 2 && split($3,arr,", ") < 2 && split($4,arr,", ") < 2 && split($5,arr,", ") < 2 && split($6,arr,", ") < 2 && split($7,arr,", ") < 2 && split($8,arr,", ") < 2 && split($9,arr,", ") < 2) print $0}' uniqorthofindOG.tsv > singlepair.tsv
# 筛选CDK基因list 不知道为什么没有筛选完成
(primer) niu 10:13:35 /mnt/c/Users/HP/Desktop/张诗诗/OG结果
$ while read line; do awk -F "\t" -v pattern="${line}" '{if(FNR==1)print $0}{if( $0 ~ pattern )print $0}' singlepair.tsv >> CDKgeneOG.tsv ; done < CDKgenelist.tsv
```
筛选出20个UTX基因，这些CDK基因在进化中有严格的同源关系。
