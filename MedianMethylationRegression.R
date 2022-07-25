library(ggplot2)

#download manifest of PMD solo-WCGWs
download.file('https://zwdzwd.s3.amazonaws.com/pmd/EPIC.comPMD.probes.tsv','./EPIC.comPMD.probes.tsv')
EPIC.comPMD<-read.delim('EPIC.comPMD.probes.tsv',header=F)

betas<-read.csv('processed.betas.csv',check.names=F,row.names=1) 
#subset betas to only PMD solo-WCGWs
b<-subset(betas,rownames(betas)%in%EPIC.comPMD$V4)
samples<-read.csv('samples.csv')
samples<-subset(samples,samples$Experiment=="Baseline profiling")
b<-b[,c(match(samples$EPIC.ID,colnames()))]
samples$med<-apply(b,2,median,na.rm=T)

cols<-c(AG21837="coral",AG06561="brown",
        AG11182="plum4",AG11546="darkseagreen4",
        AG16146="goldenrod",AG21839="darkslateblue",
        AG21859="darkslategray3")

g.PDL<-ggplot(data=samples,aes(x=Total.PDL,y=med,col=Coriell.ID))+
  geom_point(alpha=0.4)+geom_smooth(method = 'lm')+
  scale_color_manual(values=cols)+theme_classic()+
  labs(x="Population doublings in culture",y="Median PMD solo-WCGW Methylation")+
  theme(legend.title = element_blank()) +
  lims(y=c(0,1))

pdf('gg.PDL.pdf') 
g.PDL
dev.off()
