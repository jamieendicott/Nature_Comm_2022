#Compare expression of DNA methyltransferases, TETs, associated regulatory enzymes...

#import RNAseq data
library(GEOquery)
library(purrr)
library(reshape2)
library(ggplot2)


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
PCNA<-as.numeric(subset(res,rownames(res)=="PCNA"))
PCNA.norm<-apply(res,1,function(x) x/PCNA)

#vis
cols<-c(AG06561="brown",
        AG11182="plum4",AG11546="darkseagreen4",
        AG16146="goldenrod",AG21839="darkslateblue",
        AG21859="darkslategray3")
                 
dat<-p[,c(4,10)]
rownames(dat)<-p$description
dat<-cbind(dat,PCNA.norm)        
mdat<-melt(dat,id.vars = c('population_doublings','coriell_id'))
mdat$population_doublings<-as.numeric(mdat$population_doublings)
mdat$coriell_id<-as.factor(mdat$coriell_id)
                 
g<-ggplot(data=mdat,aes(x=population_doublings,y=value,col=coriell_id))
pdf('methyltransferase.exprs.pdf',width = 9,height = 3)
g+geom_point(alpha=0.5)+
  scale_color_manual(values=cols)+theme_bw()+
  labs(x="Population doublings in culture",y="log2(cpm+1)")+
  facet_wrap(~variable,ncol=7)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
                 
                 
