library(ggplot2)
library(sesame)
library(lme4)
library(multcomp)
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

#serum withdrawal
s<-subset(p,p$subexperiment=="Growth arrest: Serum withdrawal")
s$population_doublings<-as.numeric(s$population_doublings)
s$days_in_culture<-as.numeric(s$days_in_culture)
b<-b[,c(match(rownames(s),colnames(b)))]

s$med<-apply(b,2,median,na.rm=T)
levels(s$culture_pctfbs)<-c('0.50%','1%','5%','15%')

serumwd.PDL<-ggplot(data=s,aes(x=population_doublings,
                               y=med,
                               col=culture_pctfbs))+
  geom_point(alpha=0.4)+theme_classic()+geom_smooth(method='lm',alpha=0.4,colour="darkgray")+
  scale_color_manual(name="% Serum",values=c('0.50%'="deeppink4",'1%'="tomato3",'5%'="darkorange3",'15%'="goldenrod2"))+
  labs(y="Median PMD solo-WCGW methylation",x="Population Doublings")

summary(lm(s$med~s$population_doublings))
#Multiple R-squared:  0.775
#p-value: 3.081e-10

serumwd.days<-ggplot(data=s,aes(x=days_in_culture,
                                y=med,
                                col=culture_pctFBS))+
  geom_point(alpha=0.4)+theme_classic()+
  scale_color_manual(name="% Serum",values=c('0.50%'="deeppink4",'1%'="tomato3",'5%'="darkorange3",'15%'="goldenrod2"))+
  labs(y="Median PMD solo-WCGW methylation",x="Days in culture")+geom_line(stat="smooth", method="lm",alpha=0.4)+
  ylim(0.5,0.75)
  
serumwd.gc<-ggplot(data=s,aes(y=population_doublings,
                              x=days_in_culture,
                              col=culture_pctFBS))+
  geom_point(alpha=0.6)+theme_classic()+
  scale_color_manual(name="% Serum",values=c('0.50%'="deeppink4",'1%'="tomato3",'5%'="darkorange3",'15%'="goldenrod2"))+
  labs(x="Days in Culture",y="Population Doublings")+geom_line(stat="smooth", method="lm",alpha=0.4)


s$m.med<-BetaValueToMValue(s$med)
s$condition<-as.factor(s$culture_pctFBS)

res<-lmer(data=s,m.med ~ days_in_culture+condition+(1|condition),REML = FALSE)
anova(res)
summary(glht(res,mcp(condition="Tukey")))
