Sys.setenv(LOAD_ONLY='1');source('work/models.R')
z<-main[stage!='Unknown'];a<-z[,.(N=.N,Localized=sum(stage=='Localized'),Distant=sum(stage=='Distant')),by=.(sex,period,age,race)]
ref<-main[sex=='Male',.N,by=.(age,race)];ref[,w:=N/sum(N)]
res<-list();forms<-~sex*period+ns(age,knots=c(50,65,75),Boundary.knots=c(18,85))+race
for(out in c('Localized','Distant')){
 X<-model.matrix(forms,a);f<-glm.fit(X,cbind(a[[out]],a$N-a[[out]]),family=binomial());class(f)<-c('glm','lm');V<-vcov(f);stopifnot(f$converged)
 gaps<-list()
 for(pp in levels(main$period)){
  pred<-list()
  for(ss in levels(main$sex)){rr<-copy(ref);rr[,`:=`(sex=factor(ss,levels=levels(main$sex)),period=factor(pp,levels=levels(main$period)))];xx<-model.matrix(forms,rr);mu<-plogis(drop(xx%*%coef(f)));pred[[ss]]<-list(est=sum(mu*rr$w),G=colSums(xx*(rr$w*mu*(1-mu))))}
  b<-pred$Male$est-pred$Female$est;G<-pred$Male$G-pred$Female$G;se<-sqrt(drop(G%*%V%*%G));gaps[[pp]]<-list(est=b,G=G)
  res[[length(res)+1]]<-data.table(outcome=out,kind='gap',period=pp,estimate=b,lower=b-1.96*se,upper=b+1.96*se,p=2*pnorm(-abs(b/se)),male=pred$Male$est,female=pred$Female$est)
 }
 b<-gaps[[4]]$est-gaps[[1]]$est;G<-gaps[[4]]$G-gaps[[1]]$G;se<-sqrt(drop(G%*%V%*%G));res[[length(res)+1]]<-data.table(outcome=out,kind='change',period='2015-2018 vs 2000-2004',estimate=b,lower=b-1.96*se,upper=b+1.96*se,p=2*pnorm(-abs(b/se)),male=NA_real_,female=NA_real_)
}
fwrite(rbindlist(res),'outputs/tables/stage_trends.csv');print(rbindlist(res)[kind=='change'])
