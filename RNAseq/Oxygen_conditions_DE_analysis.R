#Note: RNA seq data was processed using the following pipeline: https://github.com/vari-bbc/rnaseq_workflow  

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
data_for_DE_exp$ox<-as.factor(c(rep("ambient",7,),rep("low",7)))
colnames(data_for_DE_exp)<-c("type","PDL","sample","ox")
data_for_DE_exp$PDL<-as.numeric(data_for_DE_exp$PDL)

download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197471/suppl/GSE197471_raw_cts_GEO.tsv.gz','./GSE197471_raw_cts_GEO.tsv.gz')
count_data<-read.table('GSE197471_raw_cts_GEO.tsv.gz',row.names=2)
dim(count_data)
#[1] 19384    80

count_data_exp<-count_data[,c(match(data_for_DE_exp$sample,colnames(count_data)))]

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
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
contrasts_as_str <- sapply(as.data.frame(combn(unique(paste0("ox", y$samples$ox)), 2)), function(comp) paste0(comp[1],"-",comp[2]))
contrasts <- makeContrasts(contrasts=contrasts_as_str, levels=design)
contrasts
                           
qlf <- lapply(rlang::set_names(colnames(contrasts), colnames(contrasts)), function(contrast){
  glmQLFTest(fit, contrast=contrasts[,contrast])
})
                           
invisible(lapply(names(qlf), function(contrast) {
  plotMD(qlf[[contrast]], status=decideTestsDGE(qlf[[contrast]]), values=c(1,-1), 
         col=c("red","blue"), legend="topright", hl.cex=0.6, main=contrast)
}))

res <- lapply(qlf, function(contrast) topTags(contrast, n = Inf))

lapply(res, function(contrast) table(as.data.frame(contrast)$FDR < 0.05))
#FALSE  TRUE 
#14097  1014 
       
lapply(res, function(contrast) {
  table(as.data.frame(contrast) %>% 
          dplyr::mutate(signif=FDR < 0.05, dir=ifelse(logFC>0,"up","dwn")) %>% 
          dplyr::select(signif, dir))
})
#$`oxambient-oxlow`
#       dir
#signif   dwn   up
#  FALSE 6994 7103
#  TRUE   641  373       
       
#volcano plot
volcano <- lapply(rlang::set_names(names(res),names(res)), function(contrast){
  toptable <- res[[contrast]]$table
  EnhancedVolcano::EnhancedVolcano(toptable=toptable, x="logFC", y="FDR", 
                                   lab=toptable$Uniq_syms, title=contrast, pCutoff=0.05, FCcutoff = 1, ylab="-log10(FDR)",
                                   ylim = c(0, max(-log10(toptable$FDR), na.rm = TRUE) + 0.5),
                                   xlim = c(-max(abs(toptable$logFC))-0.5, max(abs(toptable$logFC))+0.5),
    col=c('gray22', 'gray22', 'gray22', 'red3'),
    colAlpha = 0.6)
})       
       
pdf('volcano.hypoxia.nosens.pdf')
volcano
dev.off()       

#FGSEA
#Map Ensembl gene IDs to symbol. First create a mapping table.
library(fgsea)
library(DESeq2)       
counts<-count_data_exp
x.samples<-data_for_DE_exp
x.counts<-counts[,c(match(x.samples$sample,colnames(counts)))]
dds <- DESeqDataSetFromMatrix(countData = x.counts,
                              colData = x.samples,
                              design= ~ PDL + ox)
       
       
res$row <- rownames(res)
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$row, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
ens2symbol
df<-as.data.frame(res)
res2 <- inner_join(df, ens2symbol, by=c("row"="ENSEMBL"))
head(res2)
       
