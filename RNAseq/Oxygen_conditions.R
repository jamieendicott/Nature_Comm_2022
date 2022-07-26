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
library(purrr)

outdir <- "de_analysis_out_files"
dir.create(outdir)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197471')
p<-(pData(g[[1]]))

#note: this is ALL samples, not just oxygen culture condition experiment. filter:
data_for_DE_exp<-subset(p,p$'subexperiment:ch1'=="Reduced oxygen")

#remove senescence timepoints (GSM5918143, GSM5918151)
data_for_DE_exp<-subset(p,p$geo_accession!="GSM5918143" &
            p$geo_accession!="GSM5918151")

#simplify
data_for_DE_exp<-data_for_DE_exp[,c(56,62,27)]
colnames(data_for_DE_exp)<-c("ox","PDL","sample")

download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197471/suppl/GSE197471_raw_cts_GEO.tsv.gz','./GSE197471_raw_cts_GEO.tsv.gz')
count_data<-read.table('GSE197471_raw_cts_GEO.tsv.gz',row.names=2)
dim(count_data)
#[1] 19384    80

count_data_exp<-count_data[,c(match(data_for_DE_exp$description,colnames(count_data)))]

gene_names_df <- data.frame(row.names = rownames(count_data_exp))
gene_names_df$Symbol <- AnnotationDbi::mapIds(eval(as.name(orgdb)), rownames(gene_names_df), 
                                              keytype="ENSEMBL", column="SYMBOL", 
                                              multiVals="first")
gene_names_df$Uniq_syms <- scater::uniquifyFeatureNames(rownames(gene_names_df), gene_names_df$Symbol)
gene_names_df$entrez <- AnnotationDbi::mapIds(eval(as.name(orgdb)), rownames(gene_names_df), 
                                              keytype="ENSEMBL", column="ENTREZID", 
                                              multiVals="first")

testthat::expect_equal(colnames(count_data), rownames(data_for_DE))
testthat::expect_equal(rownames(count_data), rownames(gene_names_df))

y <- DGEList(count_data_exp, samples = data_for_DE_exp, genes = gene_names_df) 
design <- model.matrix(~0+PDL+ox, data = y$samples)
design

min_cpm_cutoff <- round(10/min(colSums(y$counts)/10^6), digits = 2)
min_cpm_cutoff
#[1] 0.28
min_samples <- 7
keep <- rowSums(cpm(y) > min_cpm_cutoff) >= min_samples 
table(keep)
#FALSE  TRUE 
# 4273 15111
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
norm_counts <- cpm(y, log=TRUE)







