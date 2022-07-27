#Timepoint 1: Low-PD, pre-immortalization
cell.line<-"AG06561"
s2<-subset(p,p$"subexperiment:ch1"=="Baseline profiling"&
            p$"coriell_id:ch1"==cell.line&
            p$"subculture:ch1"=="b")
dim(s2)
s2<-s2[order(s2$"population_doublings:ch1"),]
b2<-betas[,c(match(rownames(s2),colnames(betas)))]

#Timepoints 2 and 3: post-immortalization
s<-subset(samples,samples$Experiment=="TERT immortalization"&
            samples$meetsQC=='TRUE'&
            samples$Expresison.vector=='TERT')
dim(s)
s<-s[order(s$"population_doublings:ch1"),]
b<-b[,c(match(rownames(s),colnames(b)))]

s3<-rbind(s2[1,],s2[nrow(s2),],s[nrow(s),])
b3<-cbind(b2[,1],b2[,ncol(b2)],b[,ncol(b)])
rownames(b3)<-rownames(b2)
b3<-as.data.frame(b3)
colnames(b3)<-c('early.mortal','late.mortal','late.immortal')

b3$first[b3[,1]>0.7]<-"High"
b3$first[b3[,1]>0.3 & b3[,1]<0.7]<-"Intermediate"
b3$first[b3[,1]<0.3]<-"Low"
b3$first<-as.factor(b3$first)

b3$mid[b3[,2]>0.7]<-"High"
b3$mid[b3[,2]>0.3 & b3[,2]<0.7]<-"Intermediate"
b3$mid[b3[,2]<0.3]<-"Low"
b3$mid<-as.factor(b3$mid)

b3$last[b3[,3]>0.7]<-"High"
b3$last[b3[,3]>0.3 & b3[,2]<0.7]<-"Intermediate"
b3$last[b3[,3]<0.3]<-"Low"
b3$last<-as.factor(b3$last)
b3<-na.omit(b3)

dat<-b3[,4:6]
dat<-me[,c(4:6,8)]
data <- gather_set_data(dat, 1:4)
data$x <- factor(data$x, levels=c("first", "mid", "last","class"))

pdf('parsets.pdf')
ggplot(data, aes(x, id = id, split = y, value = 1)) +
  geom_parallel_sets(aes(fill = first), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'white')+
  theme_classic()
dev.off()
