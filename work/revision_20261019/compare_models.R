source('work/revision_20261019/engine.R')
ee<-new.env(parent=globalenv());expr<-parse(file.path(source_root,'work/models.R'))
for(e in expr)if(is.call(e)&&as.character(e[[1]]) %in% c('<-','=') && is.call(e[[3]]) && identical(e[[3]][[1]],as.name('function')))eval(e,ee)
form3<-~sex*period+ns(age,knots=c(50,65,75),Boundary.knots=c(18,85))+interval+sex:interval+stage+stage:interval+grade+er+race
form3flex<-~sex*period+sex*ns(age,knots=c(50,65,75),Boundary.knots=c(18,85))+interval+sex:interval+stage+stage:interval+grade+er+race
calc3<-function(z,ff,label){
 ref<-z[sex=='Male',.N,by=vars];ref[,weight:=N/sum(N)];p<-length(ff[[1]]$b);V<-matrix(0,3*p,3*p)
 for(k in 1:3){ii<-((k-1)*p+1):(k*p);V[ii,ii]<-ff[[k]]$V}
 rows<-list();add<-function(o,a,ss='Difference',age=NA_real_){se<-sqrt(pmax(0,rowSums((o$G%*%V)*o$G)));rows[[length(rows)+1]]<<-data.table(model=label,analysis=a,sex=ss,age=age,outcome=names(o$est),estimate=o$est,se=se)}
 minus<-function(m,f)list(est=m$est-f$est,G=m$G-f$G)
 m<-ee$predict_risk(ff,ref,'Male','2015-2018');f<-ee$predict_risk(ff,ref,'Female','2015-2018');add(m,'overall','Male');add(f,'overall','Female');add(minus(m,f),'overall')
 ref<-z[sex=='Male',.N,by=.(stage,grade,er,race)];ref[,`:=`(weight=N/sum(N),age=50,sex=factor('Male',levels=levels(z$sex)),period=factor('2015-2018',levels=levels(z$period)))];gs<-list()
 for(a in c(50,75)){ref[,age:=a];m<-ee$predict_risk(ff,ref,'Male','2015-2018');f<-ee$predict_risk(ff,ref,'Female','2015-2018');g<-minus(m,f);add(g,'age',age=a);gs[[as.character(a)]]<-g}
 add(minus(gs[['75']],gs[['50']]),'age75_minus50');rbindlist(rows)
}
d<-load_data();chain<-as.integer(commandArgs(TRUE)[1]);done<-character();idle<-0
while(length(done)<10){
 files<-list.files(file.path(outroot,'private_imputations'),pattern=paste0('^smc_c',chain,'_i[0-9]+[.]rds$'),full.names=TRUE)
 for(fn in setdiff(files,done)){
  label<-sub('[.]rds$','',basename(fn));dest<-file.path(outroot,'fits',paste0(label,'_model_comparison.csv'))
  if(!file.exists(dest)){
   imp<-tryCatch(readRDS(fn),error=function(e)NULL);if(is.null(imp))next
   z<-copy(d);for(v in c('stage','grade','er'))z[,(v):=imp[[v]]]
   ans<-list();for(md in c('three_shared_age','three_sex_specific_age')){
    fm<-if(md=='three_shared_age')form3 else form3flex
    ff<-ee$fit_hazards(z,fm);ans[[md]]<-calc3(z,ff,md)
   }
   fwrite(rbindlist(ans)[,imputation:=label],dest)
  }
  done<-c(done,fn);idle<-0;cat('COMPARISON',label,'\n');flush.console()
 }
 if(length(done)<10){Sys.sleep(15);idle<-idle+1;if(idle>240)stop('Timeout')}
}
