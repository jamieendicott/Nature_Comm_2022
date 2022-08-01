library(ggplot2)
library(sesame)
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

##MMC growth arrest
s<-subset(p,p$subexperiment=="Growth arrest: MMC")

#duplicating parent/NA timepoints and assigning to each condition
s$population_doublings<-as.numeric(s$population_doublings)
s$days_in_culture<-as.numeric(s$days_in_culture)
s<-s[order(s$days_in_culture),]#order by time
s2<-rbind(s[1:3,],s)
s2$culture_agent<-as.character(s2$culture_agent)
s2$culture_agent[1:3]<-"control"
s2$culture_agent[4:6]<-"MMC"

b<-b[,c(match(s2$geo_accession,colnames(b)))]
s2$med<-apply(b,2,median,na.rm=T)
s2$culture_agent<-as.factor(s2$culture_agent)

MMC.PDL<-ggplot(data=s2,aes(x=population_doublings,y=med))+geom_smooth(method='lm',alpha=0.4,colour="darkgray")+
  geom_point(aes(colour=culture_agent),alpha=0.6)+
scale_color_manual(name="Agent",values=c(control='maroon',MMC='steelblue4'))+
  theme_classic()+
  labs(y="Median PMD solo-WCGW Methylation",x="Population doublings")+
  facet_wrap(~coriell_id,scales='free_x',labeller = labeller('Cell.line' = c("AG16146"="Adult Skin Fibroblast \nAG16146",
                                                                          "AG11182" = "Adult Vascular Endothelial Cell \nAG11182",
                                                                          "AG11546" = "Adult Vascular Smooth Muscle Cell \nAG11546")))+
  theme(strip.text.x = element_text(size = 7.5))
pdf('MMC.PDL.pdf')
MMC.PDL
dev.off()

#Regression coefficients etc for each line
c<-"AG16146"
d<-subset(s2,s2$coriell_id==c)
summary(lm(d$med~d$population_doublings))
#AG11182 r2 0.7421  p-value: 7.52e-05
#AG11546 r2 0.8032 p-value: 1.437e-05
#AG16146 r2 0.9632 p-value: 5.733e-10

MMC.days<-ggplot(data=s2,aes(x=days_in_culture,y=med,col=culture_agent))+geom_point(alpha=0.6)+
   scale_color_manual(name="Agent",values=c(control='maroon',MMC='steelblue4'))+
  geom_line(stat="smooth", method="lm",alpha=0.4)+
  theme_classic()+
  labs(y="Median PMD solo-WCGW Methylation",x="Days in culture")+
  facet_wrap(~coriell_id,scales='free_x',labeller = labeller('Cell.line' = c("AG16146"="Adult Skin Fibroblast \nAG16146",
                                                                             "AG11182" = "Adult Vascular Endothelial Cell \nAG11182",
                                                                             "AG11546" = "Adult Vascular Smooth Muscle Cell \nAG11546")))+
  theme(strip.text.x = element_text(size = 7.5))



MMC.growthcurve<-ggplot(data=s2,aes(y=population_doublings,x=days_in_culture,col=culture_agent))+
  geom_point(alpha=0.6)+theme_classic()+scale_color_manual(name="Condition",values=c(control='maroon',MMC='steelblue4'))+
  labs(x="Days in culture",y="Population doublings")+geom_line(stat="smooth", method="lm",alpha=0.4)+
  facet_wrap(~coriell_id,scales='free_x',labeller = labeller('Cell.line' = c("AG16146"="Adult Skin Fibroblast \nAG16146",
                                                                            "AG11182" = "Adult Vascular Endothelial Cell \nAG11182",
                                                                            "AG11546" = "Adult Vascular Smooth Muscle Cell \nAG11546")))+
  theme(strip.text.x = element_text(size = 7.5))

s2$m.med<-BetaValueToMValue(s2$med)
s2$condition<-as.factor(s2$culture_agent)
cell.line<-"AG16146"
s<-subset(s2,s2$coriell_id==cell.line)
res<-lmer(data=s,m.med ~ days_in_culture+condition+(1|condition),REML = FALSE)
anova(res)

m<-BetaValueToMValue(b)
cell.line<-"AG11182" #etc
s3<-subset(s,s$coriell_id==cell.line)
s3 #ensure order is correct
m2<-m[,c(match(s3$geo_accession,colnames(m)))]
m.mmc<-apply(m2[,11:13],1,mean,na.rm=T)
m.c<-apply(m2[,9:10],1,mean,na.rm=T)

t.test(m2[,1],m.c,alternative='greater')
t.test(m2[,1],m.mmc,alternative='greater')
