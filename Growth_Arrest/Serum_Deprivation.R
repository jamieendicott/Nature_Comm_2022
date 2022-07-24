library(ggplot2)
library(sesame)

#serum withdrawal
betas<-read.csv('ga.sw.all.betas.csv',row.names=1,check.names = FALSE)
s<-subset(samples,samples$Experiment=="Growth arrest: Serum withdrawal")
b<-subset(betas,rownames(betas)%in%CpGs$V4)
b<-b[,c(match(s$EPIC.ID,colnames(b)))]
s$med<-apply(b,2,median,na.rm=T)
levels(s$Culture.condition..FBS)<-c('0.50%','1%','5%','15%')

serumwd.PDL<-ggplot(data=s,aes(x=Total.PDL,y=med,col=Culture.condition..FBS))+
  geom_point(alpha=0.4)+theme_classic()+geom_smooth(method='lm',alpha=0.4,colour="darkgray")+
  scale_color_manual(name="% Serum",values=c('0.50%'="deeppink4",'1%'="tomato3",'5%'="darkorange3",'15%'="goldenrod2"))+
  labs(y="Median PMD solo-WCGW methylation",x="Population Doublings")

summary(lm(s$med~s$Total.PDL))
#Multiple R-squared:  0.775
#p-value: 3.081e-10

serumwd.days<-ggplot(data=s,aes(x=Days.in.culture,y=med,col=Culture.condition..FBS))+
  geom_point(alpha=0.4)+theme_classic()+
  scale_color_manual(name="% Serum",values=c('0.50%'="deeppink4",'1%'="tomato3",'5%'="darkorange3",'15%'="goldenrod2"))+
  labs(y="Median PMD solo-WCGW methylation",x="Days in culture")+geom_line(stat="smooth", method="lm",alpha=0.4)+
  ylim(0.5,0.75)
  
serumwd.gc<-ggplot(data=s,aes(y=Total.PDL,x=Days.in.culture,col=Culture.condition..FBS))+
  geom_point(alpha=0.6)+theme_classic()+
  scale_color_manual(name="% Serum",values=c('0.50%'="deeppink4",'1%'="tomato3",'5%'="darkorange3",'15%'="goldenrod2"))+
  labs(x="Days in Culture",y="Population Doublings")+geom_line(stat="smooth", method="lm",alpha=0.4)


samples$m.med<-BetaValueToMValue(samples$med)
samples$condition<-as.factor(samples$Culture.condition..FBS)
s<-samples
res<-lmer(data=s,m.med ~ Days.in.culture+condition+(1|condition),REML = FALSE)
anova(res)
summary(glht(res,mcp(condition="Tukey")))
