#Gene expression vs PMD solo-WCGW methylation
library(GEOquery)
library(purrr)
library(ggplot2)
library(edgeR)
library(reshape2)
library(sesame)

#import manifest data
download.file('https://zhouserver.research.chop.edu/InfiniumAnnotation/20210615/EPIC/EPIC.hg38.manifest.gencode.v36.tsv.gz','./EPIC.hg38.manifest.gencode.v36.tsv.gz')
manifest.hg38 <- read.delim("EPIC.hg38.manifest.gencode.v36.tsv.gz")
download.file('https://zwdzwd.s3.amazonaws.com/pmd/EPIC.comPMD.probes.tsv','./EPIC.comPMD.probes.tsv')
CpGs<-read.delim('EPIC.comPMD.probes.tsv',header=F)
solo.man<-subset(manifest.hg38,manifest.hg38$probeID%in%CpGs$V4)

#import RNAseq data
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197471')
p<-(pData(g[[1]]))
#beautify pdata
p<-p[,c(1,2,27,55:60,62:64)]
p$population_doublings<-as.numeric(p$population_doublings)
cols<-(as.character(map(strsplit(colnames(p), split = ":"), 1)))
colnames(p)<-cols

download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197471/suppl/GSE197471_raw_cts_GEO.tsv.gz','./GSE197471_raw_cts_GEO.tsv.gz')
cts<-read.table('GSE197471_raw_cts_GEO.tsv.gz',row.names=2)
dim(cts)
#[1] 19384    80

cpm<-cpm(cts[,5:80])
cpm<-log2(cpm+1)
cpm<-cbind(cts[,1:4],cpm)

s2<-subset(p,p$subexperiment=="Baseline profiling"&
             p$coriell_id=="AG11182"|p$coriell_id=="AG21859")
e<-t(subset(cpm,cpm$Uniq_syms=="ADAMTS2"|cpm$Uniq_syms=="CARD11"))
e<-as.data.frame(e[5:80,])
colnames(e)<-c("ADAMTS2","CARD11") #can also try other genes
s2$ADAMTS2<-e[c(match(s2$description,rownames(e))),1]
s2$CARD11<-e[c(match(s2$description,rownames(e))),2]
s2$ADAMTS2<-as.numeric(s2$ADAMTS2)
s2$CARD11<-as.numeric(s2$CARD11)
s2<-s2[,c(2,4,13,14)]
s2<-melt(s2,id.vars=c('geo_accession','coriell_id'))

pdf('ADAMTS2.CARD11.log2cpm.pdf',height = 4,width = 4)
ggplot(data=s2,aes(x=variable,y=value))+geom_boxplot()+
  theme_bw()+facet_wrap(~coriell_id,scales = 'free_y')+ylim(0,10)
dev.off()

#Methylation heatmaps (SeSaMe)
#load methylation data
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

c<-subset(p,p$subexperiment=="Baseline profiling"&p$coriell_id=="AG11182") #AG21859
c<-c[order(c$population_doublings),]

b.hm<-b[,c(match(c$geo_accession,colnames(b)))]
gene<-"ADAMTS2" #CARD11, etc
visualizeGene(paste(gene), b.hm, platform = "EPIC",
                              refversion = "hg38")
