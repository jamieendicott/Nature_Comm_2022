library(purrr)
library(reshape2)
library(ggplot2)
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

##Load in methylation data
betas<-g.data
#beautify pdata...
p<-p[,c(1,2,48:56,58:61)] #can trim further if desired
cols<-(as.character(map(strsplit(colnames(p), split = ":"), 1)))
colnames(p)<-cols

cell.line<-"AG21839"
samples<-subset(p, p$subexperiment=="Baseline profiling"&
                   p$coriell_id==cell.line)
#order samples by advancing PDs
samples<-samples[order(population_doublings),]
b<-betas[,c(match(samples$geo_accession,colnames(betas)))]
dim(b)
#[1] 865918     27
delta<-b-b[,1]

names<-c('PMD.SCGS.social','HMD.SCGS.social',
           'PMD.SCGW.social','HMD.SCGW.social',
           'PMD.WCGW.social','HMD.WCGW.social',
           'PMD.SCGS.solo','HMD.SCGS.solo',
           'PMD.SCGW.solo','HMD.SCGW.solo',
           'PMD.WCGW.solo','HMD.WCGW.solo')

for (i in seq_along(names)) {
  name<-names[i]
  d<-subset(all,all$class==name)
  db<-subset(delta,rownames(delta)%in%d$probeID)
  med <- apply(db,2,median,na.rm=T)
  new_name <- paste0("med",names[i])
  assign(new_name,med)
} 

samples$PMD.SCGS.social.med<-medPMD.SCGS.social
samples$HMD.SCGS.social.med<-medHMD.SCGS.social
samples$PMD.SCGW.social.med<-medPMD.SCGW.social
samples$HMD.SCGW.social.med<-medHMD.SCGW.social
samples$PMD.WCGW.social.med<-medPMD.WCGW.social
samples$HMD.WCGW.social.med<-medHMD.WCGW.social

samples$PMD.SCGS.solo.med<-medPMD.SCGS.solo
samples$HMD.SCGS.solo.med<-medHMD.SCGS.solo
samples$PMD.SCGW.solo.med<-medPMD.SCGW.solo
samples$HMD.SCGW.solo.med<-medHMD.SCGW.solo
samples$PMD.WCGW.solo.med<-medPMD.WCGW.solo
samples$HMD.WCGW.solo.med<-medHMD.WCGW.solo

median.PDL<-samples[,c(12,16:27)]
median.PDL<-melt(median.PDL,id='population_doublings')
median.PDL$population_doublings<-as.numeric(median.PDL$population_doublings)
median.PDL$value<-as.numeric(median.PDL$value)
median.PDL$variable<-as.factor(median.PDL$variable)
median.PDL$CpG.context<-c(rep('social',162),rep('solo',162))
median.PDL$PMDvHMD<-as.factor(c(rep('PMD',27),rep('HMD',27),
                                rep('PMD',27),rep('HMD',27),
                                rep('PMD',27),rep('HMD',27),
                                rep('PMD',27),rep('HMD',27),
                                rep('PMD',27),rep('HMD',27),
                                rep('PMD',27),rep('HMD',27)))

g<-ggplot(data=median.PDL,aes(x=population_doublings,y=value,col=variable))
pdf('AG21839.noCGI.delta.PMDvHMD.pdf') 
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
