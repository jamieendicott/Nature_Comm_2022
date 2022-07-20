#Methylation array pipeline
library(sesameData)
library(sesame)
library(data.table)
#may need to cache sesame if newly installed

idat_dir = ('./idats') #path to idats
betas = do.call(cbind, lapply(searchIDATprefixes(idat_dir), function(pfx) {
  getBetas(dyeBiasNL(noob(pOOBAH(readIDATpair(pfx),pval.threshold = 0.1)))) 
})) 
#opensesame command default is p 0.05
dim(betas)

#subset PMDsoloWCGWs out
download.file('https://zwdzwd.s3.amazonaws.com/pmd/EPIC.comPMD.probes.tsv',destfile = './EPIC.comPMD.probes.tsv')
soloWCGWprobes<-read.delim('EPIC.comPMD.probes.tsv',header=F)
betas<-subset(betas,rownames(betas)%in%soloWCGWprobes$V4)
dim(betas)
