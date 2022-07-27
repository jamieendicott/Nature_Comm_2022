#Adding replication timing information to array manifest
#bigwig repliseq files from ENCODE
#note: repliseq data is mapped to hg19
library(GenomicRanges)
library(rtracklayer)
#importing repliseq .bws
#can either download locally here: http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUwRepliSeq
#or use the following
urls<-c(#BJ
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjG1bPctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjS1PctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjS2PctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjS3PctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjS4PctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjG2PctSignalRep1.bigWig',
        #HUVEC
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecG1bPctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecS1PctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecS2PctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecS3PctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecS4PctSignalRep1.bigWig',
        'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecG2PctSignalRep1.bigWig')
dests<-c('./BJ.G1b.bigWig',
         './BJ.S1.bigWig',
         './BJ.S2.bigWig',
         './BJ.S3.bigWig',
         './BJ.S4.bigWig',
         './BJ.G2.bigWig',
         './HUVEC.G1b.bigWig',
         './HUVEC.S1.bigWig',
         './HUVEC.S2.bigWig',
         './HUVEC.S3.bigWig',
         './HUVEC.S4.bigWig',
         './HUVEC.G2.bigWig')

for(i in seq_along(urls)){
    download.file(urls[i], dests[i], mode="wb")
}

RS.HUVEC.G1b<-import('HUVEC.G1b.bigWig',format="BigWig")
RS.HUVEC.S1<-import('HUVEC.S1.bigWig',format="BigWig")
RS.HUVEC.S2<-import('HUVEC.S2.bigWig',format="BigWig")
RS.HUVEC.S3<-import('HUVEC.S3.bigWig',format="BigWig")
RS.HUVEC.S4<-import('HUVEC.S4.bigWig',format="BigWig")
RS.HUVEC.G2<-import('HUVEC.G2.bigWig',format="BigWig")

RS.BJ.G1b<-import('BJ.G1b.bigWig',format="BigWig")
RS.BJ.S1<-import('BJ.S1.bigWig',format="BigWig")
RS.BJ.S2<-import('BJ.S2.bigWig',format="BigWig")
RS.BJ.S3<-import('BJ.S3.bigWig',format="BigWig")
RS.BJ.S4<-import('BJ.S4.bigWig',format="BigWig")
RS.BJ.G2<-import('BJ.G2.bigWig',format="BigWig")

#load probe manifest (PMD solo-WCGWs only here, but could use complete manifest)
EPIC.comPMD.probes <- read.delim("EPIC.comPMD.probes.tsv", header=FALSE)
solo.WCGW.universe<-GRanges(seqnames = EPIC.comPMD.probes$V1,
                            ranges = IRanges(start = EPIC.comPMD.probes$V2,
                                             end = EPIC.comPMD.probes$V3))
values(solo.WCGW.universe)<-EPIC.comPMD.probes$V4
#not all repliseq files have the same coverage, combine into overlapping GR
for(i in seq_along(grs)){
    OL.RS.HUVEC<-subsetByOverlaps(RS.HUVEC.G1b,grs[i])
}

OL.RS.HUVEC<-subsetByOverlaps(RS.HUVEC.G1b,RS.HUVEC.S1)
OL.RS.HUVEC<-subsetByOverlaps(RS.HUVEC.G1b,RS.HUVEC.S2)
OL.RS.HUVEC<-subsetByOverlaps(RS.HUVEC.G1b,RS.HUVEC.S3)
OL.RS.HUVEC<-subsetByOverlaps(RS.HUVEC.G1b,RS.HUVEC.S4)
OL.RS.HUVEC<-subsetByOverlaps(RS.HUVEC.G1b,RS.HUVEC.G2)
#...and so on
length(OL.RS.HUVEC)
#[1] 2711188
#place all RS files into universe
RS.HUVEC.G1b<-subsetByOverlaps(RS.HUVEC.G1b,OL.RS.HUVEC)
RS.HUVEC.S1<-subsetByOverlaps(RS.HUVEC.S1,OL.RS.HUVEC)
RS.HUVEC.S2<-subsetByOverlaps(RS.HUVEC.S2,OL.RS.HUVEC)
RS.HUVEC.S3<-subsetByOverlaps(RS.HUVEC.S3,OL.RS.HUVEC)
RS.HUVEC.S4<-subsetByOverlaps(RS.HUVEC.S4,OL.RS.HUVEC)
RS.HUVEC.G2<-subsetByOverlaps(RS.HUVEC.G2,OL.RS.HUVEC)
OL.RS.BJ<-subsetByOverlaps(RS.BJ.G1b,RS.BJ.S1)
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,RS.BJ.S2)
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,RS.BJ.S3)
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,RS.BJ.S4)
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,RS.BJ.G2)
length(subsetByOverlaps(OL.RS.BJ,solo.WCGW.universe))
#[1] 23762
test.cpgs<-subsetByOverlaps(test.cpgs,test.scores)
length(test.cpgs)
#[1] 25144
#lists CpGs that overlap RT file
test.scores<-subsetByOverlaps(OL.RS.BJ,solo.WCGW.universe)
length(test.scores)
#[1] 23762
#RT bins that contain CpGs, meaning there are more RT bins than CpGs
test<-findOverlapPairs(OL.RS.BJ,solo.WCGW.universe)
#where did metadata go/
test<-pintersect(findOverlapPairs(OL.RS.BJ,solo.WCGW.universe))
#cpg names left off...but scores there
test<-pintersect(findOverlapPairs(solo.WCGW.universe,OL.RS.BJ))
#now cpg names on there. combine both!
test.cpgs<-pintersect(findOverlapPairs(OL.RS.BJ.G1b.scores,solo.WCGW.universe))
test.scores<-pintersect(findOverlapPairs(solo.WCGW.universe,OL.RS.BJ))
solo.WCGW.RT<-as.data.frame(cbind(as.character(test.scores$X),test.cpgs$score))

#get all GRanges together
OL<-subsetByOverlaps(OL.RS.BJ.G1b.scores,OL.RS.BJ.S1.scores)
#...
length(OL)
#[1] 25189
#now 24559


cpgs<-pintersect(findOverlapPairs(solo.WCGW.universe,OL))

OL.RS.HUVEC.G1b.scores<-pintersect(findOverlapPairs(RS.HUVEC.G1b,solo.WCGW.universe))
OL.RS.HUVEC.S1.scores<-pintersect(findOverlapPairs(RS.HUVEC.S1,solo.WCGW.universe))
OL.RS.HUVEC.S2.scores<-pintersect(findOverlapPairs(RS.HUVEC.S2,solo.WCGW.universe))
OL.RS.HUVEC.S3.scores<-pintersect(findOverlapPairs(RS.HUVEC.S3,solo.WCGW.universe))
OL.RS.HUVEC.S4.scores<-pintersect(findOverlapPairs(RS.HUVEC.S4,solo.WCGW.universe))
OL.RS.HUVEC.G2.scores<-pintersect(findOverlapPairs(RS.HUVEC.G2,solo.WCGW.universe))

OL.RS.BJ.G1b.scores<-pintersect(findOverlapPairs(RS.BJ.G1b,solo.WCGW.universe))
OL.RS.BJ.S1.scores<-pintersect(findOverlapPairs(RS.BJ.S1,solo.WCGW.universe))
OL.RS.BJ.S2.scores<-pintersect(findOverlapPairs(RS.BJ.S2,solo.WCGW.universe))
OL.RS.BJ.S3.scores<-pintersect(findOverlapPairs(RS.BJ.S3,solo.WCGW.universe))
OL.RS.BJ.S4.scores<-pintersect(findOverlapPairs(RS.BJ.S4,solo.WCGW.universe))
OL.RS.BJ.G2.scores<-pintersect(findOverlapPairs(RS.BJ.G2,solo.WCGW.universe))

RT.soloWCGW<-cbind(as.character(cpgs$X),as.numeric(as.character(OL.RS.BJ.G1b.scores$score)))
solo.WCGW.RT$HUVEC.G1b<-OL.RS.HUVEC.G1b.scores$score


OL.RS.BJ.G1b.scores<-subsetByOverlaps(OL.RS.BJ.G1b.scores,OL)
OL.RS.BJ.S1.scores<-subsetByOverlaps(OL.RS.BJ.S1.scores,OL)
OL.RS.BJ.S2.scores<-subsetByOverlaps(OL.RS.BJ.S2.scores,OL)
OL.RS.BJ.S3.scores<-subsetByOverlaps(OL.RS.BJ.S3.scores,OL)
OL.RS.BJ.S4.scores<-subsetByOverlaps(OL.RS.BJ.S4.scores,OL)
OL.RS.BJ.G2.scores<-subsetByOverlaps(OL.RS.BJ.G2.scores,OL)
OL.RS.HUVEC.G1b.scores<-subsetByOverlaps(OL.RS.HUVEC.G1b.scores,OL)
OL.RS.HUVEC.S1.scores<-subsetByOverlaps(OL.RS.HUVEC.S1.scores,OL)
OL.RS.HUVEC.S2.scores<-subsetByOverlaps(OL.RS.HUVEC.S2.scores,OL)
OL.RS.HUVEC.S3.scores<-subsetByOverlaps(OL.RS.HUVEC.S3.scores,OL)
OL.RS.HUVEC.S4.scores<-subsetByOverlaps(OL.RS.HUVEC.S4.scores,OL)
OL.RS.HUVEC.G2.scores<-subsetByOverlaps(OL.RS.HUVEC.G2.scores,OL)

RT.soloWCGW$BJ.G1b<-as.numeric(as.character(OL.RS.BJ.G1b.scores$score))
RT.soloWCGW$BJ.S1<-as.numeric(as.character(OL.RS.BJ.S1.scores$score))
RT.soloWCGW$BJ.S2<-as.numeric(as.character(OL.RS.BJ.S2.scores$score))
RT.soloWCGW$BJ.S3<-as.numeric(as.character(OL.RS.BJ.S3.scores$score))
RT.soloWCGW$BJ.S4<-as.numeric(as.character(OL.RS.BJ.S4.scores$score))
RT.soloWCGW$BJ.G2<-as.numeric(as.character(OL.RS.BJ.G2.scores$score))

RT.soloWCGW$HUVEC.G1b<-as.numeric(as.character(OL.RS.HUVEC.G1b.scores$score))
RT.soloWCGW$HUVEC.S1<-as.numeric(as.character(OL.RS.HUVEC.S1.scores$score))
RT.soloWCGW$HUVEC.S2<-as.numeric(as.character(OL.RS.HUVEC.S2.scores$score))
RT.soloWCGW$HUVEC.S3<-as.numeric(as.character(OL.RS.HUVEC.S3.scores$score))
RT.soloWCGW$HUVEC.S4<-as.numeric(as.character(OL.RS.HUVEC.S4.scores$score))
RT.soloWCGW$HUVEC.G2<-as.numeric(as.character(OL.RS.HUVEC.G2.scores$score))



library(ggsci)
library(ggplot2)

samples<-read.csv('compiled.samples.batches.csv')
betas<-read.csv('betas.GEO.csv',header=T,row.names = 1,check.names = FALSE)

cell.line<-"AG21859"
s<-subset(samples,samples$Coriell.ID==cell.line&
                  samples$Experiment=="Baseline profiling")
s<-s[order(s$Total.PDL),]
b<-betas[,c(match(s$EPIC.ID,colnames(betas)))]

#only interested in probes that haven't already lost all methylation
b$delta<-b[,ncol(b)]-b[,1]
b<-subset(b,b[,1]>=0.3)

RTscores <- read.csv("~/Dropbox/replication timing/replication.timing/RTscores.soloWCGW.csv")
RTscores<-RTscores[!duplicated(RTscores$V2), ]
rownames(RTscores)<-RTscores$V2
RTscores<-RTscores[,-c(1,2)]
RT.rowanno<-RTscores[,c(6,5,4,3,2,1)] #BJ fibroblast
RT.rowanno<-RTscores[,c(12,11,10,9,8,7)] #HUVECs
colnames(RT.rowanno)<-c('G2','S4','S3','S2','S1','G1b')
WAscore<-as.data.frame(0.917*RT.rowanno$G1b+0.750*RT.rowanno$S1+0.583*RT.rowanno$S2+
                         0.417*RT.rowanno$S3+0.250*RT.rowanno$S4+0*RT.rowanno$G2)
WAscore[is.na(WAscore)] <- 0
colnames(WAscore)<-'score'
rownames(WAscore)<-rownames(RT.rowanno)
head(WAscore)

#regression, slope vs WAscore
test<-na.omit(b)
test<-test[,c(match(s$EPIC.ID,colnames(test)))]

lm<-apply(test,1,function(x) lm(x~s$Total.PDL)) 
fit<-lm[[1]]
names(summary(fit))
head(fit)$coefficients

B1<-vector(mode='numeric',length=length(lm))
for (i in 1:length(lm)){
  fit<-lm[[i]]
  B1[i] <- summary(fit)$coefficients["s$Total.PDL","Estimate"]
}
B1<-as.data.frame(B1)
rownames(B1)<-rownames(test)
b2<-na.omit(B1)
b2$WA<-WAscore[c(match(rownames(b2),rownames(WAscore))),1]
b2<-na.omit(b2)
cor(b2$B1,b2$WA)

#AG11182
#[1] 0.2576768
#AG21859
#[1] 0.3819171

summary(lm(b2$B1~b2$WA))
#AG11182
#Multiple R-squared:  0.0664,	Adjusted R-squared:  0.06633 
#F-statistic: 948.2 on 1 and 13332 DF,  p-value: < 2.2e-16

#AG21859
#Multiple R-squared:  0.1459,	Adjusted R-squared:  0.1458 
#F-statistic:  1553 on 1 and 9096 DF,  p-value: < 2.2e-16

#AG21859.slopevWA.all
ggplot(data=b2,aes(x=WA,y=B1))+geom_point(alpha=0.15)+theme_bw()+ggtitle('AG21859')+
  geom_smooth(method='lm',se=FALSE,col="darkred",size=0.8)

#by WA quintile
WAscore$ile<-ntile(WAscore$score,5)
WAscore$ile<-as.factor(WAscore$ile)
b2$ile<-WAscore[c(match(rownames(b2),rownames(WAscore))),2]

pdf(paste0(cell.line,'.violin.B1vWA.pdf'))
ggplot(data=b2,aes(x=ile,y=B1,fill=ile))+geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.5)+
  theme_bw()+
  ggtitle(paste0(cell.line,' slope vs WA quintile'))+
  scale_fill_manual(values = alpha(c("grey22","grey42","grey62","grey82","white"),0.2))+
  xlab('Replication timing WA score quintile')+ylab('PMD solo-WCGW methylation change per PD')
dev.off()

kruskal.test(b2$B1~b2$ile)
#Kruskal-Wallis chi-squared = 1068.6, df = 4, p-value < 2.2e-16 AG11182
#Kruskal-Wallis chi-squared = 1376.4, df = 4, p-value < 2.2e-16 AG21859
