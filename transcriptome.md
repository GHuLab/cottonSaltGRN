1. download data：
login https://sra-explorer.info/
conda activate rnaseq
cd  /vol3/agis/huguanjing_group/xiongxianpeng/salt2023/rawdata
# mkdir rawdata 

vi download.sh
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR120/012/SRR12077612/SRR12077612_1.fastq.gz . && mv SRR12077612_1.fastq.gz SRR12077612_LfJ02-508-3DPA1_1.fastq.gz
)

nohup bash download.sh &
####jobs 
###top -u xiongxianpeng
admin536363

2. #fastqc
vi fq.sh
qcdir=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/fq     

fqdir=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02   
# fqdir
/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/fastqc -t 3 -o $qcdir $fqdir/*.fastq.gz


#fastp
ls *gz |sed 's/_2.fastq.gz//g' |sed 's/_1.fastq.gz//g' | uniq>31sample.txt

vi fastp.sh

while read line
do
clean=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/clean

raw=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02

read1=${line}_1.fastq.gz
read2=${line}_2.fastq.gz
fastp -w 8 -i $raw/$read1 \
-I $raw/$read2 \
-o $clean/${line}_1.clean.gz \
-O $clean/${line}_2.clean.gz \
-l 15 -q 15 --compression=6 -R $cleandata/${line} \
-h $clean/${line}.fastp.html \
-j $clean/${line}.fastp.json
done < /vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02/31sample.txt 
qsub fastp.sh 

grep -n "total_reads"  *.json|grep 'json:15'|sed 's/.fastp.json:4://g'|sed 's/_/\t/g'>clean.read.txt
grep 'gc_content' *.json|sed '1~2d' |sed 's/"gc_content"://g'>gc.txt
grep 'q20_rate' *.json|sed '1~2d' |sed 's/"q20_rate"://g'>q20.txt
grep 'q30_rate' *.json|sed '1~2d' |sed 's/"q30_rate"://g'>q30.txt
grep 'total_bases' *.json|sed '1~2d' |awk '!(NR%2)'| sed 's/"total_bases"://g'>clean.base.txt


########Ga
#Hisat2
hisat2-build Ghirsutum_527_v2.0.fa Ghirsutum_527_v2.0 


vi hisat1.sh
while read line
do
index=/vol3/agis/huguanjing_group/xiongxianpeng/fiberData/UTXgenome/Ghirsutum_527_v2.0
inputdir=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/clean
outdir=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/Allbam
read1=${line}_1.clean.gz
read2=${line}_2.clean.gz
hisat2 -p 10 -x $index -1 $inputdir/$read1 -2 $inputdir/$read2 -S ${outdir}/${line}.Hisat_aln.sam
/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/samtools sort -@ 3 -o ${outdir}/${line}.Hisat_aln.sorted.bam ${outdir}/${line}.Hisat_aln.sam
/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/samtools view -bhF 12 -q 30 ${outdir}/${line}.Hisat_aln.sorted.bam >${outdir}/${line}.unique.bam
/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/samtools index ${outdir}/${line}.unique.bam ${outdir}/${line}.unique.bam.bai>${line}.log
/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/samtools flagstat -@ 3 ${outdir}/${line}.unique.bam > ${outdir}/${line}.txt
rm ${outdir}/${line}.Hisat_aln.sam
rm ${outdir}/${line}.Hisat_aln.sorted.bam
done < /vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02/31sample.txt 
ls *unique.bam|sed 's/.unique.bam//g'|sed 's/_/\t/g'>sample.txt
le hisat1.sh.e1494852 |grep 'reads; of these'|sed 's/ /\t/g'>read.txt
le hisat1.sh.e1494852|grep 'aligned exactly 1 time'|sed 's/ /\t/g'>unique.txt


#### 7.featurecounts
##### 7.1 featureCounts (Version 2.0.1)
vi feature.sh
inputdir=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/Allbam

gtf=/vol3/agis/huguanjing_group/xiongxianpeng/fiberData/AllBam/feature/Ghirsutum_527_v2.1.gene_exons.primaryOnly.gft3
featureCounts -T 10 -p -t CDS -g gene_id -a $gtf -o /vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/Allbam/all.id.txt $inputdir/*.unique.bam
qsub feature.sh
multiqc all.id.txt.summary
sed -i 's###g' all.id.txt
sed -i 's###g' all.id.txt
# sed 's/*_//'all.id.txt | sed -r 's/.*({7}.)/\1/' >all2.id.txt
sed -i '1d' all.id.txt
cat all.id.txt | cut -f1,7- > counts1.txt

#8.1 stringtie（4 mins per sample）
vi stringtie.sh
while read line
do
sam=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/hisat
/vol3/agis/huguanjing_group/xiongxianpeng/miniconda3/envs/rnaseq/bin/stringtie -e -B -p 8 -G /vol3/agis/huguanjing_group/xiongxianpeng/fiberData/AllBam/feature/Ghirsutum_527_v2.1.gene_exons.primaryOnly.gft3  -o /vol3/agis/huguanjing_group/xiongxianpeng/salt2023/rawdata/hisat/${line}/${line}.gtf ${sam}/${line}.unique.bam（stringtie前面要加绝对路径，which stringtie）
done </vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02/31sample.txt

#8.2  Obtain TPM
vi TPM.sh
while read line
do
bam=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/hisat
grep 'TPM' ${bam}/${line}/${line}.gtf |grep 'Goh'|cut -f 9|grep 'TPM'|awk '{print $4$10}'| sed 's/;/\t/g'|less -S> ${bam}/${line}.txt
done </vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02/31sample.txt

vi bam.sh
while read line
do
bam=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/hisat
grep 'Goh' ${bam}/${line}.txt |less -S> ${bam}/${line}-1.txt
done </vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02/31sample.txt

vi TPM2.sh
while read line
do
TPM=/vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/hisat
sort -n ${TPM}/${line}.txt > ${TPM}/${line}_TPM.sort.txt 
done < /vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02/31sample.txt

ls *sort.txt >list.txt
for i in `cat list.txt`;do cut -f2 $i >$i.tmp;done
paste -s *.tmp> all.txt
le SRR8089823_transcriptome_sequencing_of_Gossypium_hirsutum_TM-1_and_Gossypium_barbadense_Hai7124_TPM.sort.txt|cut -f 1 > a.txt
le a.txt|tr "\n" ","|sed -e 's/,$/\n/'|sed 's/,/\t/g'>aa.txt
cat  aa.txt all.txt >all.tpm.txt
paste  /vol3/agis/huguanjing_group/xiongxianpeng/PRJNA490626/haorui02/31sample.txt all.tpm.txt >finally.txt(sample +gene) 




