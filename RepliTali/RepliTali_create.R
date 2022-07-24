library('glmnet')
library('coefplot')
library('ggplot2')

cols<-c(AG21837="coral",AG06561="brown",
        AG11182="plum4",AG11546="darkseagreen4",
        AG16146="goldenrod",AG21839="darkslateblue",
        AG21859="darkslategray3")
        
man<-read.table('EPICnonCGI.context.manifest.tsv') 

samples<-read.csv('compiled.samples.batches.csv')
samples<-subset(samples,samples$Experiment=="Baseline profiling")

betas<-read.csv('Baseline profiling all.betas.csv',check.names=FALSE,row.names=1)
dim(betas)
#[1] 865918    182

b<-betas[,c(match(samples$EPIC.ID,colnames(betas)))]

#Only including PMD probes in RepliTali, subset
probeset<-"PMD"
probes<-subset(man,man$PMD=="TRUE")
dim(probes)
#56936    11
betas<-subset(b,rownames(b)%in%probes$probeID)

##Create normalization model based on chronologically youngest cell line AG06561
samples.651<-subset(samples,samples$Coriell.ID=="AG06561")
betas.651<-betas[,c(match(samples.651$EPIC.ID,colnames(betas)))]
betas.651<-na.omit(betas.651)
dim(betas.651)
y.651<-paste(samples.651$Total.PDL)
y.651<-as.numeric(y.651)
y.651
x.651<-t(betas.651)

ELR.651<-glmnet(x.651,y.651,family="gaussian",alpha=0.5) 

cv.fit.651<-cv.glmnet(x.651,y.651,alpha=0.5)
dim(extract.coef(cv.fit.651,lambda="lambda.1se")) #Note: there is a minor variability to glmnet function. size of model will vary slightly.

#Normalize starting PDLs ONLY using model built from AG06561
norm651<-extract.coef(cv.fit.651,lambda="lambda.1se")

norm651.betas<-betas[c(match(norm651$Coefficient[-1],rownames(betas))),]
norm.PDL<-apply(norm651.betas,2,function(x) norm651$Value[-1]*x)
norm.PDL<-apply(norm.PDL,2,function(x) sum(x,na.rm=T) +norm651$Value[1]) 
samples$norm.PDL<-(norm.PDL)

#for all subsequent timepoints calculate delta PDL using observed PDL
X <- split(samples, samples$Coriell.ID)
for( i in seq_along(X)){
  X[[i]]$deltaPDL<-(X[[i]]$Total.PDL-X[[i]]$Total.PDL[1])
  X[[i]]$adj651PDL<-(X[[i]]$deltaPDL+X[[i]]$norm.PDL[1])
} 
X2<-do.call("rbind", X)

##Using 651 model-adjusted starting PDL, train new model with random sample selection
set.seed(12345)
betas<-na.omit(betas) #critical step
rtrain.betas<-betas[,sample(ncol(betas),122)] #2:1 split, 60 samples in test set
rtrain.samples<-X2[c(match(colnames(rtrain.betas),X2$EPIC.ID)),]
rtest.samples<-subset(X2,!(X2$EPIC.ID%in%colnames(rtrain.betas)))
rtest.betas<-betas[,c(match(rtest.samples$EPIC.ID,colnames(betas)))]

x.r<-t(rtrain.betas)
y.r<-paste(rtrain.samples$adj651PDL)
y.r<-as.numeric(y.r)
y.r

ELR.r<-glmnet(x.r,y.r,family="gaussian",alpha=0.5)
cv.fit.ELR.r<-cv.glmnet(x.r,y.r,alpha=0.5)
dim(extract.coef(cv.fit.ELR.r,lambda="lambda.1se"))
#[1] 88   3 
cvfit.r<-extract.coef(cv.fit.ELR.r,lambda="lambda.1se")

#performance on random training set
r.enet.betas<-subset(rtrain.betas,rownames(rtrain.betas)%in%rownames(cvfit.r))
training.enet<-apply(r.enet.betas,2,function(x) cvfit.r$Value[-1]*x)
est.training<-apply(training.enet,2,function(x) sum(x) +cvfit.r$Value[1])
rtrain.samples$replitali<-paste(est.training)
#performance on random test set
r.enet.test.betas<-subset(rtest.betas,rownames(rtest.betas)%in%rownames(cvfit.r))
test.enet<-apply(r.enet.test.betas,2,function(x) cvfit.r$Value[-1]*x)
est.test<-apply(test.enet,2,function(x) sum(x) +cvfit.r$Value[1])
rtest.samples$replitali<-paste(est.test)

#combine for vis
rtrain.samples$set<-"TRAINING"
rtest.samples$set<-"TEST"
dat<-rbind(rtrain.samples,rtest.samples)
dat$replitali<-as.numeric(dat$replitali)

g<-ggplot(data=transform(dat,set=factor(set,levels=c("TRAINING","TEST"))),aes(x=adj651PDL,y=mitoage))
g+geom_point(aes(col=Coriell.ID),size=1.5,alpha=0.6)+theme_bw()+
  scale_color_manual(values=cols)+
  geom_smooth(method='lm', alpha = 0.4,col='gray')+
  xlab("Observed PDs, 651 enet adjusted P1") + ylab("Predicted PDs")+facet_wrap(~set)

