#Heatmap of AG21859 PMD solo-WCGWs, showing reproducibility between replicates/subcultures 
library(pheatmap)
library(viridisLite)

#download manifest of PMD solo-WCGWs
download.file('https://zwdzwd.s3.amazonaws.com/pmd/EPIC.comPMD.probes.tsv','./EPIC.comPMD.probes.tsv')
EPIC.comPMD<-read.delim('EPIC.comPMD.probes.tsv',header=F)

betas<-read.csv('processed.betas.csv',check.names=F,row.names=1) 
#subset betas to only PMD solo-WCGWs
b<-subset(betas,rownames(betas)%in%EPIC.comPMD$V4)

s<-subset(samples,samples$Coriell.ID=="AG21859")
#sort by subculture then PDs
s<-s[order( s[,4], s[,5] ),]
#move first timepoint to first row
s2<-rbind(s[40,],s[-40,])

b<-betas[,c(match(s2$EPIC.ID,colnames(betas)))]
b<-na.omit(b)
b2<-b[order(b[,1],decreasing = TRUE),]

p<-pheatmap(b3,cluster_cols = F, cluster_rows =F,
         show_rownames = F, show_colnames = F,
         gaps_col = c(14,26),
         color=turbo(100)
)

pdf('hm.AG21859.subs.pdf')
p
dev.off()
