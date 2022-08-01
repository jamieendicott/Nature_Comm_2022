#load methylation data
library(GEOquery)
library(purrr)
library(ggplot2)
library(ggforce)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197512')
p<-pData(g[[1]])
betas <- as.data.frame(exprs(g[[1]]))
dim(betas)
#[1] 865918    372
#beautify pdata
p<-p[,c(1,2,48:56,58:60)]
cols<-(as.character(map(strsplit(colnames(p), split = ":"), 1)))
colnames(p)<-cols
#subset betas to PMD solo-WCGWs
download.file('https://zwdzwd.s3.amazonaws.com/pmd/EPIC.comPMD.probes.tsv','./EPIC.comPMD.probes.tsv')
EPIC.comPMD<-read.delim('EPIC.comPMD.probes.tsv',header=F)
b<-subset(betas,rownames(betas)%in%EPIC.comPMD$V4)
dim(b)
#[1] 26732   372

#Timepoint 1: Low-PD, pre-immortalization
#Timepoint 2: High-PD, pre-immortalization
cell.line<-"AG06561"
s2<-subset(p,p$subexperiment=="Baseline profiling"&
            p$coriell_id==cell.line&
            p$subculture=="b")
dim(s2)
s2<-s2[order(s2$population_doublings),]
b2<-b[,c(match(rownames(s2),colnames(b)))]

#Timepoint 3: post-immortalization
s<-subset(p,p$subexperiment=="TERT immortalization"&
            p$expression_vector=='TERT')
dim(s)
s<-s[order(s$population_doublings),]
b<-b[,c(match(rownames(s),colnames(b)))]

s3<-rbind(s2[1,],s2[nrow(s2),],s[nrow(s),])
b3<-cbind(b2[,1],b2[,ncol(b2)],b[,ncol(b)])
rownames(b3)<-rownames(b2)
b3<-as.data.frame(b3)
colnames(b3)<-c('early.mortal','late.mortal','late.immortal')

b3$first[b3[,1]>0.7]<-"High"
b3$first[b3[,1]>0.3 & b3[,1]<0.7]<-"Intermediate"
b3$first[b3[,1]<0.3]<-"Low"
b3$first<-as.factor(b3$first)

b3$mid[b3[,2]>0.7]<-"High"
b3$mid[b3[,2]>0.3 & b3[,2]<0.7]<-"Intermediate"
b3$mid[b3[,2]<0.3]<-"Low"
b3$mid<-as.factor(b3$mid)

b3$last[b3[,3]>0.7]<-"High"
b3$last[b3[,3]>0.3 & b3[,2]<0.7]<-"Intermediate"
b3$last[b3[,3]<0.3]<-"Low"
b3$last<-as.factor(b3$last)
b3<-na.omit(b3)

dat<-b3[,4:6]
data <- gather_set_data(dat, 1:3) 
data$x <- factor(data$x, levels=c("first", "mid", "last","class"))

pdf('parsets.pdf')
ggplot(data, aes(x, id = id, split = y, value = 1)) +
  geom_parallel_sets(aes(fill = first), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'white')+
  theme_classic()
dev.off()
