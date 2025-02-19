library(ggplot2)
c <- read.table(file = "C:/Users/hr345/Desktop/all/new/Öù×´Í¼.txt", sep = "\t", header = T )
head(c)
#group1 group2 number
#1     CK     up   2151
#2     CK   down    596
#3      S     up    560
#4      S   down     93
#5      A     up    241
#6      A   down    253
ggplot( c, aes( x = group1, y = number, fill = group2))+ geom_bar(stat="identity", color="black")+scale_fill_brewer(palette="Pastel")


##YYµÄÖù×´Í¼
rm(list=ls())
library(ggplot2)
library(dplyr)
library(magrittr)
c <- read.table(file = "C:/Users/hr345/Desktop/´æ»îÂÊ.txt", sep = "\t", header = T )
rate_mean <- c %>% 
     dplyr::group_by(variety) %>% 
     dplyr::summarize(
         count=n(),
         mean = mean(rate),
         sd = sd(rate))
plot_data1 <- rate_mean
plot_data2 <- c
p4 <- ggplot()+ 
     geom_bar(data=plot_data1,mapping=aes(x=variety,y=mean,fill=variety),position="dodge",stat="identity",width = 0.7)+geom_errorbar(data=plot_data1,mapping=aes(x = variety,ymin = mean-sd, ymax = mean+sd),width = 0.1,color = 'black', size=0.8)+theme_classic(base_line_size = 1 )+labs(y="Survival Rate(%)")+scale_fill_manual(values=c("#5e5e5f","#3f4040","#272625","#444f59","#3e4149","#1e2022"))+ theme(legend.position="none")
p4         