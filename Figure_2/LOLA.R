library(transcripTools)
library(dplyr)
library(ggplot2)
library("LOLA")
library('GenomicRanges')
library(simpleCache)

s<-subset(p,p$"subexperiment:ch1"=="TERT immortalization"&
            p$"expression_vector:ch1"=="TERT")
dim(s)
s<-s[order(s$"population_doublings:ch1"),]

b<-betas[,c(match(rownames(s),colnames(betas)))]
dim(b)
#[1] 26732    65
b$delta<-b[,ncol(b)]-b[,1]

hi.me<-subset(b,b[,1]>0.7) 
hi.me<-subset(hi.me,abs(hi.me$delta)<0.1) 
dim(hi.me)
lo.me<-subset(b,b[,1]<0.3)
lo.me<-subset(lo.me,abs(lo.me$delta)<0.1)
var.me<-subset(b,abs(b$delta)>0.1)
#split variable into quartiles based on starting methylation
var.me$ile<-ntile(var.me[,1],4)
v1<-subset(var.me,var.me$ile==1)
v2<-subset(var.me,var.me$ile==2)
v3<-subset(var.me,var.me$ile==3)
v4<-subset(var.me,var.me$ile==4)

U<-read.delim('soloWCGW.universe.bed',header=F)
regionDB<-loadRegionDB('regions/LOLACore/hg19')
setwd('/secondary/projects/laird/jamie/LOLA/AG06561/')
CpGs<-read.table('EPIC.comPMD.probes.tsv')
#convert to granges object, can also be granges list
userUniverse<-GRanges(seqnames = U$V1,
                      ranges = IRanges(start = U$V2,
                                       end = U$V3))

group<-v4 #etc
input<-subset(CpGs,CpGs$V4%in%rownames(group))
input<-input[,1:3]
userSets<-GRanges(seqnames = input$V1,
                  ranges = IRanges(start = input$V2,
                                   end = input$V3))
locResults = runLOLA(userSets, userUniverse, regionDB, cores=1)
head(locResults)

writeCombinedEnrichment(locResults, outFolder= "lolaRes.TERT.v4")
E<-read.delim('allEnrichments.v4.tsv')#etc
#select top 10
E<-E[1:10,c(2,4,5,15)]
#annotate description w dbset
E$Set<-paste(E$description," db",E$dbSet,sep = "")
E$group<-'low'
Elo<-E

E<-rbind(Ehi,EV4,EV3,EV2,EV1,Elo)
E$group<-as.factor(E$group)
E$Set<-reorder(E$Set, E$pValueLog,mean)

g<-ggplot(data=E,aes(x=pValueLog,y=Set,fill=oddsRatio))
pdf('TERT.enrichment.LOLA.pdf',width = 8,height=8)
g+geom_col()+theme_bw()+
  xlab('-log(p-val)')+
  scale_fill_distiller(palette = "Spectral")+
  facet_grid(rows = vars(group),scales = "free_y")
dev.off()
