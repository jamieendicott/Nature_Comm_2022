man<-read.table('EPICnonCGI.context.manifest.tsv') 

PMD.SCGS.social<-subset(man,man$PMD==TRUE & man$SCGS==TRUE & man$type=='social') #7418
PMD.SCGS.social$class<-'PMD.SCGS.social'
HMD.SCGS.social<-subset(man,man$HMD==TRUE & man$SCGS==TRUE & man$type=='social') #4691
HMD.SCGS.social$class<-'HMD.SCGS.social'
PMD.SCGW.social<-subset(man,man$PMD==TRUE & man$SCGW==TRUE & man$type=='social') #4566
PMD.SCGW.social$class<-'PMD.SCGW.social'
HMD.SCGW.social<-subset(man,man$HMD==TRUE & man$SCGW==TRUE & man$type=='social') #2900
HMD.SCGW.social$class<-'HMD.SCGW.social'
PMD.WCGW.social<-subset(man,man$PMD==TRUE & man$WCGW==TRUE & man$type=='social') #2804
PMD.WCGW.social$class<-'PMD.WCGW.social'
HMD.WCGW.social<-subset(man,man$HMD==TRUE & man$WCGW==TRUE & man$type=='social') #1804
HMD.WCGW.social$class<-'HMD.WCGW.social'

PMD.SCGS.solo<-subset(man,man$PMD==TRUE & man$SCGS==TRUE & man$type=='solo') #4470
PMD.SCGS.solo$class<-'PMD.SCGS.solo'
HMD.SCGS.solo<-subset(man,man$HMD==TRUE & man$SCGS==TRUE & man$type=='solo') #3103
HMD.SCGS.solo$class<-'HMD.SCGS.solo'
PMD.SCGW.solo<-subset(man,man$PMD==TRUE & man$SCGW==TRUE & man$type=='solo') #5939
PMD.SCGW.solo$class<-'PMD.SCGW.solo'
HMD.SCGW.solo<-subset(man,man$HMD==TRUE & man$SCGW==TRUE & man$type=='solo') #3796
HMD.SCGW.solo$class<-'HMD.SCGW.solo'
PMD.WCGW.solo<-subset(man,man$PMD==TRUE & man$WCGW==TRUE & man$type=='solo') #8356
PMD.WCGW.solo$class<-'PMD.WCGW.solo'
HMD.WCGW.solo<-subset(man,man$HMD==TRUE & man$WCGW==TRUE & man$type=='solo') #5068
HMD.WCGW.solo$class<-'HMD.WCGW.solo'

all<-rbind(PMD.SCGS.social,HMD.SCGS.social,
           PMD.SCGW.social,HMD.SCGW.social,
           PMD.WCGW.social,HMD.WCGW.social,
           PMD.SCGS.solo,HMD.SCGS.solo,
           PMD.SCGW.solo,HMD.SCGW.solo,
           PMD.WCGW.solo,HMD.WCGW.solo)
           
for (i in names) {
  name<-names[i]
  d<-subset(all,all$class==name)
  db<-subset(delta,rownames(delta)%in%d$probeID)
  med <- apply(db,2,median,na.rm=T)
  assign(paste0(names[i],'.med'),med)
} 

name<-'HMD.WCGW.solo'
d<-subset(all,all$class==name)
db<-subset(betas,rownames(betas)%in%d$probeID)
samples$HMD.WCGW.solo<-apply(db,2,median)

median.PDL<-median.betas[,c(2,4:15)]
median.PDL<-melt(median.PDL,id='Total.PDL')
median.PDL$Total.PDL<-as.numeric(median.PDL$Total.PDL)
median.PDL$value<-as.numeric(median.PDL$value)
median.PDL$variable<-as.factor(median.PDL$variable)
c<-median.PDL
c$CpG.context<-c(rep('social',162),rep('solo',162))
c$PMDvHMD<-c(rep('PMD',27),rep('HMD',27),rep('PMD',27),rep('HMD',27),rep('PMD',27),rep('HMD',27),rep('PMD',27),rep('HMD',27),rep('PMD',27),rep('HMD',27),rep('PMD',27),rep('HMD',27))

g<-ggplot(data=c,aes(x=Total.PDL,y=value,col=variable))
pdf('AG21839.noCGI.delta.PMDvHMD.pdf') #also made one with methylation
g+geom_point()+stat_smooth(geom='line',alpha=0.4,method='lm')+
  theme_classic()+xlab('Total PDL in Culture')+
  facet_wrap(~PMDvHMD)+
  ylab('Delta median methylation')+
  scale_color_manual(values = c('firebrick','firebrick',
                                'aquamarine4','aquamarine4',
                                'skyblue4','skyblue4',
                                'firebrick1','firebrick1',
                                'aquamarine2','aquamarine2',
                                'skyblue1','skyblue1'))
  
dev.off()
