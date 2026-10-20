source('work/revision_20261019/engine.R')
d<-load_data();files<-list.files(file.path(outroot,'private_imputations'),pattern='^smc_.*[.]rds$',full.names=TRUE)
stopifnot(length(files)==40);cv<-c('stage','grade','er','race');ans<-list()
for(fn in files){
 z<-copy(d);imp<-readRDS(fn);for(v in c('stage','grade','er'))z[,(v):=imp[[v]]]
 ref<-z[sex=='Male',.N,by=cv];ref[,weight:=N/sum(N)]
 counts<-z[,.N,by=c('sex','period','age',cv)]
 for(scope in c('Full fitted cohort','2015-2018 only'))for(a in 40:80){
  s<-counts[age>=a-5 & age<a+5 & (scope=='Full fitted cohort' | period=='2015-2018')]
  ct<-dcast(s[,.(N=sum(N)),by=c('sex',cv)],stage+grade+er+race~sex,value.var='N',fill=0)
  ct<-merge(ref,ct,by=cv,all.x=TRUE);ct[is.na(Male),Male:=0];ct[is.na(Female),Female:=0]
  ans[[length(ans)+1]]<-data.table(imputation=basename(fn),scope=scope,age=a,zero_weight=sum(ct[Male==0 | Female==0]$weight),below5_weight=sum(ct[Male<5 | Female<5]$weight),male_n=sum(s[sex=='Male']$N),female_n=sum(s[sex=='Female']$N))
 }
}
x<-rbindlist(ans);fwrite(x,file.path(outroot,'support_grid_individual.csv'))
fwrite(x[,.(zero_weight=mean(zero_weight),below5_weight=mean(below5_weight),male_n=mean(male_n),female_n=mean(female_n)),by=.(scope,age)],file.path(outroot,'support_grid.csv'))
