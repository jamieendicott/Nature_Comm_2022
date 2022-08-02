download.file('https://raw.githubusercontent.com/jamieendicott/Nature_Comm_2022/main/RepliTali/RepliTali_coefs.csv','./RepliTali_coefs.csv')
RT<-read.csv('RepliTali_coefs.csv')
#start with betas object, sort to match order with samples of interest in separate df
RT.betas<-subset(betas,rownames(betas)%in%RT$Coefficient)
res-apply(RT.betas,2,function(x) RT$Value[-1]*x)
res<-apply(res,2,function(x) sum(x) +RT$Value[1])
