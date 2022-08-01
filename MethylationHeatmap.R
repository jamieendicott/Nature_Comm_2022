#Heatmap of AG21859 PMD solo-WCGWs, showing reproducibility between replicates/subcultures 
library(pheatmap)
library(viridisLite)
#load methylation data
library(GEOquery)
library(purrr)
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
s<-subset(p,p$subexperiment=="Baseline profiling" &
          p$coriell_id=="AG21859")
#sort by subculture then PDs
s<-s[order( s[,13], s[,12] ),]


b<-b[,c(match(rownames(s),colnames(b)))]
b<-na.omit(b)
b2<-b[order(b[,1],decreasing = TRUE),]

p<-pheatmap(b2,cluster_cols = F, cluster_rows =F,
         show_rownames = F, show_colnames = F,
         gaps_col = c(14,26),
         color=turbo(100)
)

pdf('hm.AG21859.subs.pdf')
p
dev.off()
