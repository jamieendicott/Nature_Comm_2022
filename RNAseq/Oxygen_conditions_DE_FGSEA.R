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

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197471')
p<-(pData(g[[1]]))

#note: this is ALL samples, not just oxygen culture condition experiment. filter:
p<-subset(p,p$'subexperiment:ch1'=="Reduced oxygen")

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
library(apeglm)
library(msigdbr)
       
       
counts<-count_data_exp
       
# set organism for loking up gene sets
msigdb_organism <- "Homo sapiens"

# set seed for gsea
gsea_seed <- 12345

de_objs <- list(dgelist=y, fit=fit, qlf=qlf, res=res)
#de_objs <- readRDS("../template_Rmd_for_downstream_work/de_analysis_out_files_oxygen/edgeR.rds")
de_res <- lapply(de_objs$res, as.data.frame)
sapply(de_res, nrow)

cpms<- gene_names_df %>% 
  rownames_to_column("ensembl_id") %>% 
  left_join(., as_tibble(norm_counts, rownames="ensembl_id"), by="ensembl_id")        
       
       
# Remove genes with 0 logFC and PValue = 1, calculate ranking metric then sort Entrez genes in descending order
# Adapts code from
# https://github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList.
prep_clusterprofiler_genelist <- function(dge_table, rank_by="signed-log10pval"){
  
  # filter away genes with exactly 0 logFC and PValue = 1.
  filt_dge_table <- dge_table[dge_table$logFC != 0 &
                                dge_table$PValue != 1, ]
  
  # calculate rank_metric
  if(identical(rank_by, "signed-log10pval")){
    filt_dge_table$rank_metric <-
      sign(filt_dge_table$logFC) * -log10(filt_dge_table$PValue)
  } else{
    stop("Specify valid ranking metric.")
  }
  
  ## feature 1: numeric vector
  geneList <- filt_dge_table$rank_metric
  ## feature 2: named vector
  names(geneList) <- as.character(filt_dge_table$entrez)
  ## feature 3: decreasing order
  geneList <- sort(geneList, decreasing = TRUE)

  return(geneList)
}


# Get the genesets and format for clusterprofiler (dataframe with col1 = geneset name, col2 = entrez gene)
# organisms is 'Homo sapiens' or 'Mus musculus'
# if no msigdb subcat, then specify NA
get_geneset <- function(gene_set, msigdb_subcat=NA, organism){
  if (gene_set %in% c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")){
    #browser()
    msigdbr_args <- list(species = organism, category = gene_set, subcat=msigdb_subcat)
    msigdbr_args <-  msigdbr_args[!sapply(msigdbr_args, is.na)] # remove 'subcat' param if it is NA
    
    msigdbr_gene_set <- do.call(msigdbr::msigdbr, msigdbr_args)
    
    # convert to clusterprofiler friendly format
    geneset_out <- msigdbr_gene_set[, c("gs_name", "entrez_gene")] %>%
      as.data.frame(stringsAsFactors = FALSE)
    
  } else{
    stop("Invalid value for gene_set parameter.")
  }
  
  geneset_out
}

# Get ranked genes
genes_and_score <- lapply(de_res, prep_clusterprofiler_genelist)
sapply(genes_and_score, length)

# remove genes with no Entrez ID (isNA) 
genes_and_score <- lapply(genes_and_score, function(x) x[!is.na(names(x))])
sapply(genes_and_score, length)
#drops down to 13765

# remove genes with duplicated Entrez ID
genes_and_score <- lapply(genes_and_score, function(x) {
  ids <- names(x)
  duplicated_ids <- unique(ids[duplicated(ids)])
  x[!ids %in% duplicated_ids]
})
sapply(genes_and_score, length)

# Confirm genes are ordered in decreasing order
correct_ranking <- sapply(genes_and_score, function(x) {
  all(order(x, decreasing = TRUE) == 1:length(x))
})
stopifnot(all(correct_ranking))

# Get genesets and run GSEA
genesets_of_interest <- list(H=c("H",NA))
genesets <- lapply(genesets_of_interest, function(x) get_geneset(gene_set=x[1], 
                                                                 msigdb_subcat=x[2], 
                                                                 organism=msigdb_organism))

gsea_res <- lapply(genesets, function(geneset){
  
  lapply(genes_and_score, function(y){
    set.seed(gsea_seed) # make reproducible
    
    # gene_list is named vector where names are the Entrez IDs and values are the ranking metric
    gsea_res <- clusterProfiler::GSEA(geneList = y,
                                      TERM2GENE = geneset, 
                                      eps = 0.0 # need to set this or Pvalues will not reach below 1e-10
                                      )
    
    gsea_res_syms <- DOSE::setReadable(gsea_res,
                                       OrgDb = eval(as.symbol(orgdb)),
                                       keyType = "ENTREZID")
    list(entrez=gsea_res, symbols=gsea_res_syms)
    
  })
})

# output to file
invisible(lapply(names(gsea_res), function(geneset){
  lapply(names(gsea_res[[geneset]]), function(contrast){
    gseaResult <- gsea_res[[geneset]][[contrast]]$symbols
    write_tsv(as_tibble(gseaResult), paste0(outdir, "/", contrast, "_", geneset, ".tsv")) 
  })
}))

#Summary plot
lapply(names(gsea_res), function(geneset){
  lapply(names(gsea_res[[geneset]]), function(contrast){
    gseaResult <- gsea_res[[geneset]][[contrast]]$symbols
    if(nrow(gseaResult) > 0){
      dotplot(gseaResult, split=".sign") + ggtitle(paste0(contrast," -- ",geneset)) + 
        scale_y_discrete(label=function(x) str_trunc(x, 40)) + facet_grid(.~.sign)
    } else{
      "No significant results" 
    }
  })
})
