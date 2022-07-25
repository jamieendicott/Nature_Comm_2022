library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE197512')
p<-pData(g[[1]])
p<-read.csv('samples.csv',row.names=1)

getGEOSuppFiles(GEO='GSE197512', makeDirectory = TRUE, baseDir = getwd(),
  fetch_files = TRUE) 
list.files(".")

