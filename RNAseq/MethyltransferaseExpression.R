#Compare expression of DNA methyltransferases, TETs, associated regulatory enzymes...

#import RNAseq data
library(GEOquery)
library(purrr)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197471')
p<-(pData(g[[1]]))
#beautify pdata
p<-p[,c(1,2,27,55:60,62:64)]
cols<-(as.character(map(strsplit(colnames(p), split = ":"), 1)))
colnames(p)<-cols
p<-subset(p,p$subexperiment=="Baseline profiling")

download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197471/suppl/GSE197471_norm_cts_GEO.tsv.gz','./GSE197471_norm_cts_GEO.tsv.gz')
cts<-read.table('GSE197471_norm_cts_GEO.tsv.gz',row.names=2)
dim(cts)
#[1] 19384    80
ncts<-cts[,c(match(p$description,colnames(cts)))]  
ncts<-cbind(cts[,1:4],ncts)

#Select genes
set<-c("DNMT1","DNMT3A","DNMT3B","UHRF1","TET1","TET2","TET3","PCNA")
res<-subset(ncts,ncts$Symbol %in% set)
rownames(res)<-res$Symbol
res<-res[,-(1:4)]

#normalize to PCNA
tres<-t(res)
PCNA<-as.numeric(subset(res,rownames(res)=="PCNA"))

PCNA.norm<-apply(tres,1,function(x) x/PCNA)
PCNA.norm<-apply(res,2,function(x) x/PCNA)
