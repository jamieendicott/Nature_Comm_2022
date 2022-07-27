library(reshape2)
library(viridisLite)
library(pheatmap)
library(ggplot2)
library(dplyr)

#subset betas to PMD solo-WCGWs
cell.line<-"AG06561"
s<-subset(p,p$"subexperiment:ch1"=="Baseline profiling"&
            p$"coriell_id:ch1"==cell.line&
            p$"subculture:ch1"=="b")
dim(s)
s<-s[order(s$"population_doublings:ch1"),]

b<-betas[,c(match(rownames(s),colnames(betas)))]
b$delta<-b[,ncol(b)]-b[,1]

hi.me<-subset(b,b[,1]>0.7) 
hi.me<-subset(hi.me,abs(hi.me$delta)<0.1) 
lo.me<-subset(b,b[,1]<0.3)
lo.me<-subset(lo.me,abs(lo.me$delta)<0.1)
var.me<-subset(b,abs(b$delta)>0.1)

#make order for heatmap
o.hime<-hi.me[order(-hi.me[,1]),]
o.lome<-lo.me[order(lo.me[,1]),]
o.varme<-var.me[order(-var.me[,1]),]
hm.order<-c(rownames(o.hime),rownames(o.varme),rownames(o.lome))

an.group<-c(rep("Low stable",length(rownames(lo.me))),
                    rep("Variable",length(rownames(var.me))),
                    rep("High stable",length(rownames(hi.me))))
an.group<-as.data.frame(an.group)
row.names(an.group)<- c(rownames(lo.me),rownames(var.me),rownames(hi.me))
an.group$an.group<-as.factor(an.group$an.group)
colnames(an.group)<-"Group"

dat<-b[,-ncol(b)] #remove delta col
dat<-na.omit(dat)
an.group<-subset(an.group,rownames(an.group)%in%rownames(dat))
dat<-subset(dat, rownames(dat)%in%rownames(an.group))
hm.order<-subset(hm.order,hm.order%in%rownames(an.group))

png(paste0(cell.line,"hm.groups.png"),width = 7,height = 5,units = 'in', res = 600)
pheatmap(dat[c(match(hm.order,rownames(dat))),],cluster_cols = F, cluster_rows =F,
         show_rownames = F, show_colnames = F,
         annotation_row = an.group, color=turbo(100))
dev.off()

##Aurora/viridis density plots
#Stably methylated/unmethylated groups
dat<-as.data.frame(hi.me)
dat<-dat[,-ncol(dat)]
colnames(dat)<-s$"population_doublings:ch1"
dat$CpG<-rownames(dat)
mdat<-melt(dat,id.vars = 'CpG') 
mdat$variable<-as.numeric(as.character(mdat$variable))
mdat$value<-as.numeric(mdat$value)
mdat<-na.omit(mdat)
head(mdat)

base_plot <- ggplot(data=mdat,aes(x=variable,y=value)) + geom_point(alpha=0)
base_plot+stat_density_2d(
  aes(fill = after_stat(..density..)),
  geom = "raster",
  contour = FALSE)+
  theme_classic()+scale_fill_viridis_c(option="A")+
  ylim(0,1)+
  xlab("PDs")+ylab("Methylation")


#variable group, further split into quartile
dat<-as.data.frame(var.me)
dat<-dat[,-ncol(dat)]
colnames(dat)<-s$"population_doublings:ch1"
dat$CpG<-rownames(dat)
#determine decile for starting methylation
dat$ile<-ntile(dat[,1],4) 
mdat<-melt(dat,id.vars = c('CpG','ile')) 
mdat$variable<-as.numeric(as.character(mdat$variable))
mdat$value<-as.numeric(mdat$value)
mdat<-na.omit(mdat)
head(mdat)

base_plot <- ggplot(data=mdat,aes(x=variable,y=value)) + geom_point(alpha=0)
base_plot+stat_density_2d(
  aes(fill = after_stat(..density..)),
  geom = "raster",
  contour = FALSE)+
  theme_classic()+scale_fill_viridis_c(option="A")+
  ylim(0,1)+
  xlab("PDs")+ylab("Methylation")+
  ggtitle("Methylation by starting value quartile")+
  facet_wrap(~ile)
