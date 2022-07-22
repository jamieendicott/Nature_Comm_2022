#profiling cpg contexts
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)  

#download complete manifest
download.file('https://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg19.manifest.tsv.gz','./man.hg19.tsv.gz')
hg19.man<-read.delim('man.hg19.tsv.gz',stringsAsFactors=FALSE)
##removing all CGI probes
EPIC.nonCGI<-subset(hg19.man,is.na(gene)) #159497
##simplify
EPIC.nonCGI<-EPIC.nonCGI[,c(1:3,5,15)] 
dim(EPIC.nonCGI)
#[1] 159497      5

g<-GRanges(seqnames = EPIC.nonCGI$CpG_chrm,
  ranges = IRanges(start = EPIC.nonCGI$CpG_beg,
    end = EPIC.nonCGI$CpG_end,
    names = EPIC.nonCGI$probeID))
g$names <- EPIC.nonCGI$probeID

##PMD vs HMD annotation
download.file('https://zwdzwd.s3.amazonaws.com/pmd/PMD_coordinates_hg19.bed.gz','./PMD_coordinates_hg19.bed.gz')
PMD_coordinates_hg19.bed <- read.delim("PMD_coordinates_hg19.bed.gz", header=FALSE, stringsAsFactors=FALSE)
PMD<-subset(PMD_coordinates_hg19.bed,PMD_coordinates_hg19.bed$V6=='commonPMD') #13127 PMDs
PMD.g<-GRanges(seqnames = PMD$V1,
  ranges=IRanges(start=PMD$V2,
  end=PMD$V3))
HMD<-subset(PMD_coordinates_hg19.bed,PMD_coordinates_hg19.bed$V6=='commonHMD') #7089 HMDs
HMD.g<-GRanges(seqnames = HMD$V1,
  ranges=IRanges(start=HMD$V2,
  end=HMD$V3))
  
PMDss<-subsetByOverlaps(g,PMD.g) #56936 CpGs
HMDss<-subsetByOverlaps(g,HMD.g) #37187 CpGs

PMDss<-as.data.frame(PMDss)
EPIC.nonCGI$PMD<-EPIC.nonCGI$probeID %in% PMDss$names
HMDss<-as.data.frame(HMDss)
EPIC.nonCGI$HMD<-EPIC.nonCGI$probeID %in% HMDss$names


##Immediate sequence context
#3 sequence contexts
#SCGS SCGW WCGW (S=C|G, W=T|A)
#load all CpGs in hg19
chrs <- names(Hsapiens)[1:24]

#SCGS context
CCGC <- lapply(chrs, function(x) start(matchPattern("CCGC", Hsapiens[[x]])))
CCGC <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(CCGC[[x]], width = 4))))
length(CCGC) #2055755
GCGG <- lapply(chrs, function(x) start(matchPattern("GCGG", Hsapiens[[x]])))
GCGG <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(GCGG[[x]], width = 4))))
length(GCGG) #2054015
GCGC <- lapply(chrs, function(x) start(matchPattern("GCGC", Hsapiens[[x]])))
GCGC <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(GCGC[[x]], width = 4))))
length(GCGC) #1655892
CCGG <- lapply(chrs, function(x) start(matchPattern("CCGG", Hsapiens[[x]])))
CCGG <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(CCGG[[x]], width = 4))))
length(CCGG) #2297198
SCGS<-c(CCGC,CCGG,GCGC,GCGG) #8062860
SCGSss<-subsetByOverlaps(g,SCGS) #44142 non-CGI SCGS on EPIC array
SCGSss<-as.data.frame(SCGSss)
EPIC.nonCGI$SCGS<-EPIC.nonCGI$probeID %in% SCGSss$names

#SCGW context
CCGA <- lapply(chrs, function(x) start(matchPattern("CCGA", Hsapiens[[x]])))
CCGA <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(CCGA[[x]], width = 4))))
length(CCGA) #1730239
CCGT <- lapply(chrs, function(x) start(matchPattern("CCGT", Hsapiens[[x]])))
CCGT <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(CCGT[[x]], width = 4))))
length(CCGT) #1817336
GCGA <- lapply(chrs, function(x) start(matchPattern("GCGA", Hsapiens[[x]])))
GCGA <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(GCGA[[x]], width = 4))))
length(GCGA) #1473831
GCGT <- lapply(chrs, function(x) start(matchPattern("GCGT", Hsapiens[[x]])))
GCGT <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(GCGT[[x]], width = 4))))
length(GCGT) #1629769
SCGW<-c(CCGA,CCGT,GCGA,GCGT) #6651175
SCGWss<-subsetByOverlaps(g,SCGW) #39203 non-CGI SCGS on EPIC array
SCGWss<-as.data.frame(SCGWss)
EPIC.nonCGI$SCGW<-EPIC.nonCGI$probeID %in% SCGWss$names

#WCGW context
ACGA <- lapply(chrs, function(x) start(matchPattern("ACGA", Hsapiens[[x]])))
ACGA <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(ACGA[[x]], width = 4))))
length(ACGA) #1592419
ACGT <- lapply(chrs, function(x) start(matchPattern("ACGT", Hsapiens[[x]])))
ACGT <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(ACGT[[x]], width = 4))))
length(ACGT) #2151075
TCGA <- lapply(chrs, function(x) start(matchPattern("TCGA", Hsapiens[[x]])))
TCGA <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(TCGA[[x]], width = 4))))
length(TCGA) #1512838
TCGT <- lapply(chrs, function(x) start(matchPattern("TCGT", Hsapiens[[x]])))
TCGT <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(TCGT[[x]], width = 4))))
length(TCGT) #1601215
WCGW<-c(ACGA,ACGT,TCGA,TCGT) #6857547
WCGWss<-subsetByOverlaps(g,WCGW) #41048 non-CGI SCGS on EPIC array
WCGWss<-as.data.frame(WCGWss)
EPIC.nonCGI$WCGW<-EPIC.nonCGI$probeID %in% WCGWss$names

##Solo CpGs vs Social CpGs
EPIC.nonCGI$type[EPIC.nonCGI$context35>=3]<-'social' #50719
EPIC.nonCGI$type[EPIC.nonCGI$context35==0]<-'solo' #68126

write.table(EPIC.nonCGI,'EPICnonCGI.context.manifest.tsv')
