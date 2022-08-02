#Important: reprocess beta values using SeSaMe pipeline
#produce object: betas (df)
library(GEOquery)
library(purrr)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
g<-getGEO('GSE179847')
p<-pData(g[[1]])
#Beautiful dataset, however for simplicity only interested in control culture conditions/WT primary cells and corresponding contact inhibition
p<-p[,c(1,2,8,51,68,70:73,76,79,83,86,88,91,97,98)] #simplify
cols<-(as.character(map(strsplit(colnames(p), split = ":"), 1)))
colnames(p)<-cols
s<-subset(p,p$clinical_condition=="Normal" &
          p$treatments=="Control" | p$treatments=="Contact_Inhibition")
s<-subset(s,s$percent_oxygen=="21")
dim(s)
#[1] 120  17
b<-betas[,c(match(rownames(s),colnames(betas)))]
