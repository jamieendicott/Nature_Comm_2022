#Note: RNA seq data was processed using the following pipeline: 

library(GEOquery)
orgdb <- "org.Hs.eg.db"
library(orgdb, character.only=TRUE) # load the org.db for your organism
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(edgeR)
library(tibble)
library(kableExtra)
library(readr)
library(stringr)
library(cowplot)
library(EnhancedVolcano)
library(patchwork)

g<-getGEO('GSE197471')
p<-(pData(g[[1]]))

download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197471/suppl/GSE197471_norm_cts_GEO.tsv.gz','./GSE197471_norm_cts_GEO.tsv.gz')
norm_cts<-read.table('GSE197471_norm_cts_GEO.tsv.gz')
dim(norm_cts)
#[1] 19384    80
#note: this is ALL samples, not just oxygen culture condition experiment. filter:
p<-subset(p,p$'subexperiment:ch1'=="Reduced oxygen")

#remove senescence timepoints (GSM5918143, GSM5918151)
p<-subset(p,p$geo_accession!="GSM5918143" &
            p$geo_accession!="GSM5918151")
            
count_data_exp<-cbind(norm_cts[,1:4],norm_cts[,c(match(p$description,colnames(norm_cts)))])

