#Adding replication timing information to array manifest
#bigwig repliseq files from ENCODE
#note: repliseq data is mapped to hg19
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
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
         #HUVEC
         './HUVEC.G1b.bigWig',
         './HUVEC.S1.bigWig',
         './HUVEC.S2.bigWig',
         './HUVEC.S3.bigWig',
         './HUVEC.S4.bigWig',
         './HUVEC.G2.bigWig')

for(i in seq_along(urls)){
    download.file(urls[i], dests[i], mode="wb")
}

RS.BJ.G1b<-import('BJ.G1b.bigWig',format="BigWig")
RS.BJ.S1<-import('BJ.S1.bigWig',format="BigWig")
RS.BJ.S2<-import('BJ.S2.bigWig',format="BigWig")
RS.BJ.S3<-import('BJ.S3.bigWig',format="BigWig")
RS.BJ.S4<-import('BJ.S4.bigWig',format="BigWig")
RS.BJ.G2<-import('BJ.G2.bigWig',format="BigWig")

RS.HUVEC.G1b<-import('HUVEC.G1b.bigWig',format="BigWig")
RS.HUVEC.S1<-import('HUVEC.S1.bigWig',format="BigWig")
RS.HUVEC.S2<-import('HUVEC.S2.bigWig',format="BigWig")
RS.HUVEC.S3<-import('HUVEC.S3.bigWig',format="BigWig")
RS.HUVEC.S4<-import('HUVEC.S4.bigWig',format="BigWig")
RS.HUVEC.G2<-import('HUVEC.G2.bigWig',format="BigWig")

#load probe manifest (PMD solo-WCGWs only here, but could use complete manifest)
download.file('https://zwdzwd.s3.amazonaws.com/pmd/EPIC.comPMD.probes.tsv','./EPIC.comPMD.probes.tsv')
EPIC.comPMD.probes <- read.delim("EPIC.comPMD.probes.tsv", header=FALSE)
solo.WCGW.universe<-GRanges(seqnames = EPIC.comPMD.probes$V1,
                            ranges = IRanges(start = EPIC.comPMD.probes$V2,
                                             end = EPIC.comPMD.probes$V3))
values(solo.WCGW.universe)<-EPIC.comPMD.probes$V4
#not all repliseq files have the same coverage, combine into overlapping GR
OL.RS.HUVEC<-subsetByOverlaps(RS.HUVEC.G1b,RS.HUVEC.S1)
OL.RS.HUVEC<-subsetByOverlaps(OL.RS.HUVEC,RS.HUVEC.S2)
OL.RS.HUVEC<-subsetByOverlaps(OL.RS.HUVEC,RS.HUVEC.S3)
OL.RS.HUVEC<-subsetByOverlaps(OL.RS.HUVEC,RS.HUVEC.S4)
OL.RS.HUVEC<-subsetByOverlaps(OL.RS.HUVEC,RS.HUVEC.G2)
length(OL.RS.HUVEC)
#[1] 2711188
#further filtered down to PMD solo-WCGWs
OL.RS.HUVEC<-subsetByOverlaps(OL.RS.HUVEC,solo.WCGW.universe)
#24543 PMD solo-WCGW with coverage on all RT files
#place all RS files into trimmed, shared set
RS.HUVEC.G1b<-subsetByOverlaps(RS.HUVEC.G1b,OL.RS.HUVEC)
RS.HUVEC.S1<-subsetByOverlaps(RS.HUVEC.S1,OL.RS.HUVEC)
RS.HUVEC.S2<-subsetByOverlaps(RS.HUVEC.S2,OL.RS.HUVEC)
RS.HUVEC.S3<-subsetByOverlaps(RS.HUVEC.S3,OL.RS.HUVEC)
RS.HUVEC.S4<-subsetByOverlaps(RS.HUVEC.S4,OL.RS.HUVEC)
RS.HUVEC.G2<-subsetByOverlaps(RS.HUVEC.G2,OL.RS.HUVEC)

OL.RS.HUVEC.G1b.scores<-pintersect(findOverlapPairs(RS.HUVEC.G1b,solo.WCGW.universe))
OL.RS.HUVEC.S1.scores<-pintersect(findOverlapPairs(RS.HUVEC.S1,solo.WCGW.universe))
OL.RS.HUVEC.S2.scores<-pintersect(findOverlapPairs(RS.HUVEC.S2,solo.WCGW.universe))
OL.RS.HUVEC.S3.scores<-pintersect(findOverlapPairs(RS.HUVEC.S3,solo.WCGW.universe))
OL.RS.HUVEC.S4.scores<-pintersect(findOverlapPairs(RS.HUVEC.S4,solo.WCGW.universe))
OL.RS.HUVEC.G2.scores<-pintersect(findOverlapPairs(RS.HUVEC.G2,solo.WCGW.universe))

cpgs<-pintersect(findOverlapPairs(solo.WCGW.universe,OL.RS.HUVEC))
#tidy up into df
probeID<-cpgs$X
G1b<-OL.RS.HUVEC.G1b.scores$score
S1<-OL.RS.HUVEC.S1.scores$score
S2<-OL.RS.HUVEC.S2.scores$score
S3<-OL.RS.HUVEC.S3.scores$score
S4<-OL.RS.HUVEC.S4.scores$score
G2<-OL.RS.HUVEC.G2.scores$score
HUVEC.RT.soloWCGW<-data.frame(probeID,G1b,S1,S2,S3,S4,G2)

#same process for BJ
OL.RS.BJ<-subsetByOverlaps(RS.BJ.G1b,RS.BJ.S1)
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,RS.BJ.S2)
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,RS.BJ.S3)
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,RS.BJ.S4)
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,RS.BJ.G2)
length(OL.RS.BJ)
#[1] 2616253
OL.RS.BJ<-subsetByOverlaps(OL.RS.BJ,solo.WCGW.universe)
#place all RS files into trimmed, shared set
RS.BJ.G1b<-subsetByOverlaps(RS.BJ.G1b,OL.RS.BJ)
RS.BJ.S1<-subsetByOverlaps(RS.BJ.S1,OL.RS.BJ)
RS.BJ.S2<-subsetByOverlaps(RS.BJ.S2,OL.RS.BJ)
RS.BJ.S3<-subsetByOverlaps(RS.BJ.S3,OL.RS.BJ)
RS.BJ.S4<-subsetByOverlaps(RS.BJ.S4,OL.RS.BJ)
RS.BJ.G2<-subsetByOverlaps(RS.BJ.G2,OL.RS.BJ)

OL.RS.BJ.G1b.scores<-pintersect(findOverlapPairs(RS.BJ.G1b,solo.WCGW.universe))
OL.RS.BJ.S1.scores<-pintersect(findOverlapPairs(RS.BJ.S1,solo.WCGW.universe))
OL.RS.BJ.S2.scores<-pintersect(findOverlapPairs(RS.BJ.S2,solo.WCGW.universe))
OL.RS.BJ.S3.scores<-pintersect(findOverlapPairs(RS.BJ.S3,solo.WCGW.universe))
OL.RS.BJ.S4.scores<-pintersect(findOverlapPairs(RS.BJ.S4,solo.WCGW.universe))
OL.RS.BJ.G2.scores<-pintersect(findOverlapPairs(RS.BJ.G2,solo.WCGW.universe))

cpgs<-pintersect(findOverlapPairs(solo.WCGW.universe,OL.RS.BJ))
probeID<-cpgs$X
G1b<-OL.RS.BJ.G1b.scores$score
S1<-OL.RS.BJ.S1.scores$score
S2<-OL.RS.BJ.S2.scores$score
S3<-OL.RS.BJ.S3.scores$score
S4<-OL.RS.BJ.S4.scores$score
G2<-OL.RS.BJ.G2.scores$score
BJ.RT.soloWCGW<-data.frame(probeID,G1b,S1,S2,S3,S4,G2)

#load in sample info, betas if haven't already
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197512')
p<-pData(g[[1]])
betas <- as.data.frame(exprs(g[[1]]))
#beautify pdata
p<-p[,c(1,2,48:56,58:60)] #can trim further if desired
cols<-(as.character(map(strsplit(colnames(p), split = ":"), 1)))
colnames(p)<-cols

##Obtain regression coefficients by CpG, by cell type
cell.line<-"AG11182" #AG21859
samples<-subset(p, p$subexperiment=="Baseline profiling"&
                  p$coriell_id==cell.line)
#order samples by advancing PDs
samples<-samples[order(samples$population_doublings),]
b<-betas[,c(match(samples$geo_accession,colnames(betas)))]
dim(b)
#[1] 865918     40 | 31

RT.rowanno<-BJ.RT.soloWCGW #BJ fibroblast compare to AG21859 fibroblast 
RT.rowanno<-HUVEC.RT.soloWCGW #HUVECs compare to AG11182 endothelial cells
RT.rowanno<-RT.rowanno[!duplicated(RT.rowanno$probeID), ]
WAscore<-as.data.frame(0.917*RT.rowanno$G1b+0.750*RT.rowanno$S1+0.583*RT.rowanno$S2+
                         0.417*RT.rowanno$S3+0.250*RT.rowanno$S4+0*RT.rowanno$G2)
colnames(WAscore)<-'score'
rownames(WAscore)<-RT.rowanno$probeID
head(WAscore)

#Regress methylation across PMD solo-WCGWs to population doublings
b<-na.omit(subset(b,rownames(b)%in%EPIC.comPMD.probes$V4))
samples$population_doublings<-as.numeric(samples$population_doublings)
lm<-apply(b,1,function(x) lm(x~samples$population_doublings)) 
fit<-lm[[1]]
names(summary(fit))
head(fit)$coefficients

B1<-vector(mode='numeric',length=length(lm))
for (i in 1:length(lm)){
  fit<-lm[[i]]
  B1[i] <- summary(fit)$coefficients["samples$population_doublings","Estimate"]
}
B1<-as.data.frame(B1)
rownames(B1)<-rownames(b)
B1$WA<-WAscore[c(match(rownames(B1),rownames(WAscore))),1]
B1<-na.omit(B1)

#Classify CpG methylation change by which WA quintile CpG is in
B1$ile<-as.factor(ntile(B1$WA,5))

pdf(paste0(cell.line,'.violin.B1vWA.pdf'))
ggplot(data=B1,aes(x=ile,y=B1,fill=ile))+geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.5)+
  theme_bw()+
  ggtitle(paste0(cell.line,' slope vs WA quintile'))+
  scale_fill_manual(values = alpha(c("grey22","grey42","grey62","grey82","white"),0.2))+
  xlab('Replication timing WA score quintile')+ylab('PMD solo-WCGW methylation change per PD')
dev.off()

kruskal.test(B1$B1~B1$ile)
#Kruskal-Wallis chi-squared = 845.03, df = 4, p-value < 2.2e-16 AG11182
#Kruskal-Wallis chi-squared = 1301.2, df = 4, p-value < 2.2e-16 AG21859
