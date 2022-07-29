library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197512')
p<-pData(g[[1]])
g.data <- as.data.frame(exprs(g[[1]]))
head(g.data)
dim(g.data)
#[1] 865918    372
