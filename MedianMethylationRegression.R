library(ggplot2)
#load methylation data
library(GEOquery)
library(purrr)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197512')
p<-pData(g[[1]])
betas <- as.data.frame(exprs(g[[1]]))
dim(betas)
#[1] 865918    372
#subset betas to PMD solo-WCGWs
download.file('https://zwdzwd.s3.amazonaws.com/pmd/EPIC.comPMD.probes.tsv','./EPIC.comPMD.probes.tsv')
EPIC.comPMD<-read.delim('EPIC.comPMD.probes.tsv',header=F)
b<-subset(betas,rownames(betas)%in%EPIC.comPMD$V4)
dim(b)
#[1] 26732   372

#beautify pdata
p<-p[,c(1,2,48:56,58:60)]
cols<-(as.character(map(strsplit(colnames(p), split = ":"), 1)))
colnames(p)<-cols

samples<-subset(p,p$subexperiment=="Baseline profiling")
b<-b[,c(match(rownames(samples),colnames(b)))]
samples$med<-apply(b,2,median,na.rm=T)

cols<-c(AG21837="coral",AG06561="brown",
        AG11182="plum4",AG11546="darkseagreen4",
        AG16146="goldenrod",AG21839="darkslateblue",
        AG21859="darkslategray3")

g.PDL<-ggplot(data=samples,aes(x=population_doublings,y=med,col=coriell_id))+
  geom_point(alpha=0.4)+geom_smooth(method = 'lm')+
  scale_color_manual(values=cols)+theme_classic()+
  labs(x="Population doublings in culture",y="Median PMD solo-WCGW Methylation")+
  theme(legend.title = element_blank()) +
  lims(y=c(0,1))

pdf('gg.PDL.pdf') 
g.PDL
dev.off()
