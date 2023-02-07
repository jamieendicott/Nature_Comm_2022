#Methylation array pipeline
#updated 02072023 as some functions have depreciated
library(sesameData)
library(sesame)
library(data.table)
#need to cache sesame if newly installed

idat_dir = ('./idats') #path to idats
betas = openSesame(idat.dir)
#opensesame command default is p 0.05, tuning this to 0.1 is recommended but stringent cutoff is acceptable here
dim(betas)

#subset PMDsoloWCGWs out
download.file('https://zhouserver.research.chop.edu/InfiniumAnnotation/EPIC/EPIC.comPMD.probes.tsv.gz',
              './EPIC.comPMD.probes.tsv')
soloWCGWprobes<-read.delim('EPIC.comPMD.probes.tsv',header=F)
betas<-subset(betas,rownames(betas)%in%soloWCGWprobes$V4)
dim(betas)
