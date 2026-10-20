source('work/revision_20261019/engine.R')
args<-commandArgs(TRUE);mode<-if(length(args))args[1] else 'baseline'
dir.create(file.path(outroot,'fits'),showWarnings=FALSE)
d<-load_data()
analyze_one<-function(z,label,intersection=TRUE){
 fn<-file.path(outroot,'fits',paste0(label,'.rds'));if(file.exists(fn))return()
 tm<-proc.time()[3];ff<-fit_models(z)
 rr<-results(ff,z,label)
 # Restricted reference supported at both ages in both sexes, with >=5 cases
 # per profile in each 10-year window; the same restricted reference is used twice.
 cv<-c('stage','grade','er','race');ref<-z[sex=='Male',.N,by=cv];ref[,weight:=N/sum(N)]
 support<-list();sets<-list();latest_sets<-list()
 for(a in c(50,75)){
  ct<-dcast(z[age>=a-5 & age<a+5,.N,by=c('sex',cv)],stage+grade+er+race~sex,value.var='N',fill=0)
  ct<-merge(ref,ct,by=cv,all.x=TRUE);ct[is.na(Male),Male:=0];ct[is.na(Female),Female:=0]
  sets[[as.character(a)]]<-ct[Male>=5 & Female>=5,cv,with=FALSE]
  lc<-dcast(z[period=='2015-2018' & age>=a-5 & age<a+5,.N,by=c('sex',cv)],stage+grade+er+race~sex,value.var='N',fill=0)
  lc<-merge(ref,lc,by=cv,all.x=TRUE);lc[is.na(Male),Male:=0];lc[is.na(Female),Female:=0]
  latest_sets[[as.character(a)]]<-lc[Male>=5 & Female>=5,cv,with=FALSE]
  support[[length(support)+1]]<-data.table(model=label,age=a,reference_weight_any_missing=sum(ct[Male==0 | Female==0]$weight),reference_weight_below5=sum(ct[Male<5 | Female<5]$weight),latest_zero_weight=sum(lc[Male==0 | Female==0]$weight),latest_below5_weight=sum(lc[Male<5 | Female<5]$weight),male_events=sum(z$sex=='Male' & z$age>=a-5 & z$age<a+5 & z$time<=5 & z$detail %in% causes))
 }
 common<-merge(sets[['50']],sets[['75']],by=cv)
 ref<-merge(ref,common,by=cv);retained<-sum(ref$weight);ref[,weight:=weight/sum(weight)]
 ref[,`:=`(age=50,sex=factor('Male',levels=levels(z$sex)),period=factor('2015-2018',levels=levels(z$period)))]
 if(nrow(ref)){
  V<-make_V(ff);gg<-list()
  for(a in c(50,75)){
   ref[,age:=a];m<-predict_joint(ff,ref,'Male');f<-predict_joint(ff,ref,'Female');g<-list(est=m$est-f$est,G=m$G-f$G);gg[[as.character(a)]]<-g
   se<-sqrt(pmax(0,rowSums((g$G%*%V)*g$G)))
   rr<-rbind(rr,data.table(model=label,analysis='supported_reference',age=a,sex='Difference',period='2015-2018',outcome=names(g$est),estimate=g$est,se=se))
  }
  G<-gg[['75']]$G-gg[['50']]$G;ee<-gg[['75']]$est-gg[['50']]$est
  rr<-rbind(rr,data.table(model=label,analysis='supported_age75_minus50',age=NA_real_,sex='Difference',period='2015-2018',outcome=names(ee),estimate=ee,se=sqrt(pmax(0,rowSums((G%*%V)*G)))))
 }
 latest_common<-merge(latest_sets[['50']],latest_sets[['75']],by=cv)
 lr<-merge(z[sex=='Male',.N,by=cv],latest_common,by=cv);latest_retained<-sum(lr$N)/sum(z$sex=='Male')
 if(nrow(lr)){
  lr[,`:=`(weight=N/sum(N),age=50,sex=factor('Male',levels=levels(z$sex)),period=factor('2015-2018',levels=levels(z$period)))]
  gs<-list();V<-make_V(ff)
  for(a in c(50,75)){
   lr[,age:=a];m<-predict_joint(ff,lr,'Male');f<-predict_joint(ff,lr,'Female');gs[[as.character(a)]]<-list(est=m$est-f$est,G=m$G-f$G)
  }
  G<-gs[['75']]$G-gs[['50']]$G;ee<-gs[['75']]$est-gs[['50']]$est
  rr<-rbind(rr,data.table(model=label,analysis='latest_supported_age75_minus50',age=NA_real_,sex='Difference',period='2015-2018',outcome=names(ee),estimate=ee,se=sqrt(pmax(0,rowSums((G%*%V)*G)))))
 }
 # Numerical derivative check on the joint integration including cross-cause effects.
 ref[,age:=60];V<-make_V(ff);o<-predict_joint(ff,ref,'Male');p<-length(ff[[1]]$b);err<-numeric()
 for(j in c(1,10,p+10,6*p+1)){
  k<-(j-1)%/%p+1;ix<-(j-1)%%p+1;fp<-fm<-ff;fp[[k]]$b[ix]<-fp[[k]]$b[ix]+1e-5;fm[[k]]$b[ix]<-fm[[k]]$b[ix]-1e-5
  fd<-(predict_joint(fp,ref,'Male')$est-predict_joint(fm,ref,'Male')$est)/2e-5
  err<-c(err,max(abs(fd-o$G[,j])))
 }
 stopifnot(max(err)<1e-6)
 # Intersect ER-positive with localized/regional, rather than separate subgroups.
 ixcounts<-NULL
 if(intersection){
  zz<-droplevels(z[er=='Positive' & stage %in% c('Localized','Regional')]);ixcounts<-zz[,.N,by=sex]
  # ER is constant in this intersection: omit it from the analysis design.
  oldform<-form;form<<-update(form,~.-er)
  fs<-fit_models(zz);rr<-rbind(rr,results(fs,zz,paste0(label,'_ERpos_LR'),curves=FALSE)[,analysis:='intersection'])
  form<<-oldform
 }
 saveRDS(list(results=rr,support=rbindlist(support),supported_reference_weight=retained,latest_supported_reference_weight=latest_retained,intersection_counts=ixcounts,gradient_error=max(err),fits=ff,n=nrow(z)),fn)
 fwrite(rr,file.path(outroot,'fits',paste0(label,'.csv')))
 cat(format(Sys.time()),'DONE',label,'seconds',proc.time()[3]-tm,'supported_weight',retained,'\n');flush.console()
}
if(mode=='baseline'){
 analyze_one(d,'unknown_category')
 analyze_one(d[stage!='Unknown' & grade!='Unknown' & er!='Unknown'],'complete_case')
 # Descriptive unadjusted Aalen-Johansen estimates, distinct from standardization.
 suppressPackageStartupMessages(library(survival));obs<-list()
 for(sx in c('Male','Female')){
  zz<-d[period=='2015-2018' & sex==sx];status<-factor(ifelse(zz$time<=5 & zz$detail %in% causes,zz$detail,'censored'),levels=c('censored',causes))
  sf<-survfit(Surv(pmin(zz$time,5),status)~1);ii<-max(which(sf$time<=5))
  for(k in seq_along(sf$states))obs[[length(obs)+1]]<-data.table(sex=sx,outcome=sf$states[k],estimate=sf$pstate[ii,k],se=sf$std.err[ii,k],n=nrow(zz))
 }
 fwrite(rbindlist(obs),file.path(outroot,'observed_AJ.csv'))
 fwrite(d[time<=5 & detail %in% causes,.N,by=.(sex,period,age,detail)],file.path(outroot,'event_support.csv'))
}
if(mode=='mi'){
 chain<-as.integer(args[2]);done<-character();idle<-0
 while(length(done)<10){
  files<-list.files(file.path(outroot,'private_imputations'),pattern=paste0('^smc_c',chain,'_i[0-9]+[.]rds$'),full.names=TRUE)
  for(fn in setdiff(files,done)){
   imp<-tryCatch(readRDS(fn),error=function(e)NULL);if(is.null(imp))next
   z<-copy(d);stopifnot(identical(z$row,imp$row));for(v in c('stage','grade','er'))z[,(v):=imp[[v]]]
   analyze_one(z,sub('[.]rds$','',basename(fn)));done<-c(done,fn);idle<-0
  }
  if(length(done)<10){Sys.sleep(15);idle<-idle+1;if(idle>240)stop('No imputations received for one hour')}
 }
}
