library(ggplot2)
library(sesame)

CpGs<-read.table('EPIC.comPMD.probes.tsv')
samples<-read.csv('compiled.samples.batches.csv')
betas<-read.csv('ga.mmc.all.betas.csv',row.names=1,check.names = FALSE)

##MMC growth arrest
s<-subset(samples,samples$Experiment=="Growth arrest: MMC")

#duplicating parent/NA timepoints and assigning to each condition
s2<-rbind(s[1:3,],s)
s2$Agent<-as.character(s2$Agent)
s2$Agent[1:3]<-"control"
s2$Agent[4:6]<-"MMC"

b<-subset(betas,rownames(betas)%in%CpGs$V4)
b<-b[,c(match(s2$EPIC.ID,colnames(b)))]
s2$med<-apply(b,2,median,na.rm=T)
s2$Agent<-as.factor(s2$Agent)

MMC.PDL<-ggplot(data=s2,aes(x=Total.PDL,y=med))+geom_smooth(method='lm',alpha=0.4,colour="darkgray")+
  geom_point(aes(colour=Agent),alpha=0.6)+scale_color_manual(name="Agent",values=c(control='maroon',MMC='steelblue4'))+
  theme_classic()+
  labs(y="Median PMD solo-WCGW Methylation",x="Population doublings")+
  facet_wrap(~Coriell.ID,scales='free_x',labeller = labeller('Cell.line' = c("AG16146"="Adult Skin Fibroblast \nAG16146",
                                                                          "AG11182" = "Adult Vascular Endothelial Cell \nAG11182",
                                                                          "AG11546" = "Adult Vascular Smooth Muscle Cell \nAG11546")))+
  theme(strip.text.x = element_text(size = 7.5))
pdf('MMC.PDL.pdf')
MMC.PDL
dev.off()

#Regression coefficients etc for each line
c<-"AG16146"
d<-subset(samples,samples$Coriell.ID==c)
summary(lm(d$med~d$Total.PDL))
#AG11182 r2 0.7421  p-value: 7.52e-05
#AG11546 r2 0.8032 p-value: 1.437e-05
#AG16146 r2 0.9632 p-value: 5.733e-10

MMC.days<-ggplot(data=s2,aes(x=Days.in.culture,y=med,col=Agent))+geom_point(alpha=0.6)+
   scale_color_manual(name="Agent",values=c(control='maroon',MMC='steelblue4'))+
  geom_line(stat="smooth", method="lm",alpha=0.4)+
  theme_classic()+
  labs(y="Median PMD solo-WCGW Methylation",x="Days in culture")+
  facet_wrap(~Coriell.ID,scales='free_x',labeller = labeller('Cell.line' = c("AG16146"="Adult Skin Fibroblast \nAG16146",
                                                                             "AG11182" = "Adult Vascular Endothelial Cell \nAG11182",
                                                                             "AG11546" = "Adult Vascular Smooth Muscle Cell \nAG11546")))+
  theme(strip.text.x = element_text(size = 7.5))



MMC.growthcurve<-ggplot(data=s2,aes(y=Total.PDL,x=Days.in.culture,col=Agent))+
  geom_point(alpha=0.6)+theme_classic()+scale_color_manual(name="Condition",values=c(control='maroon',MMC='steelblue4'))+
  labs(x="Days in culture",y="Population doublings")+geom_line(stat="smooth", method="lm",alpha=0.4)+
  facet_wrap(~Coriell.ID,scales='free_x',labeller = labeller('Cell.line' = c("AG16146"="Adult Skin Fibroblast \nAG16146",
                                                                            "AG11182" = "Adult Vascular Endothelial Cell \nAG11182",
                                                                            "AG11546" = "Adult Vascular Smooth Muscle Cell \nAG11546")))+
  theme(strip.text.x = element_text(size = 7.5))

samples$m.med<-BetaValueToMValue(samples$med)
samples$condition<-as.factor(samples$Agent)
cell.line<-"AG16146"
s<-subset(samples,samples$Coriell.ID==cell.line)
res<-lmer(data=s,m.med ~ Days.in.culture+condition+(1|condition),REML = FALSE)
anova(res)

s3<-subset(s2,s2$Coriell.ID==cell.line)
s3 #ensure order is correct
m2<-m[,c(match(s3$EPIC.ID,colnames(m)))]
m.mmc<-apply(m2[,12:14],1,mean,na.rm=T)
m.c<-apply(m2[,10:11],1,mean,na.rm=T)

t.test(m2[,1],m.c,alternative='greater'))
t.test(m2[,1],m.mmc,alternative='greater'))


