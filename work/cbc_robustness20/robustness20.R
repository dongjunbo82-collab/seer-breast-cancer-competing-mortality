source('work/revision_20261019/engine.R')
dest<-'outputs/cbc_robustness20';dir.create(dest,recursive=TRUE,showWarnings=FALSE)
args<-commandArgs(TRUE);chain<-as.integer(args[1]);limit<-if(length(args)>1)as.integer(args[2]) else 10L
d<-load_data()
specs<-list(two_knots=c(55,70),shifted_knots=c(45,60,75))
for(ii in seq_len(limit)){
 label<-sprintf('smc_c%d_i%02d',chain,ii);imp<-readRDS(file.path(outroot,'private_imputations',paste0(label,'.rds')))
 stopifnot(identical(d$row,imp$row));z<-copy(d);for(v in c('stage','grade','er'))z[,(v):=imp[[v]]]
 ref<-z[sex=='Male',.N,by=.(stage,grade,er,race)];ref[,`:=`(weight=N/sum(N),age=50,sex=factor('Male',levels=levels(z$sex)),period=factor('2015-2018',levels=levels(z$period)))]
 for(spec in names(specs)){
  fn<-file.path(dest,paste0(spec,'_',label,'.csv'));if(file.exists(fn))next
  tm<-proc.time()[3];kn<-paste(specs[[spec]],collapse=',')
  form<-as.formula(paste0('~sex*period+sex*ns(age,knots=c(',kn,'),Boundary.knots=c(18,85))+I((age-65)*(sex=="Male"))+interval+sex:interval+stage+stage:interval+grade+er+race'))
  ff<-fit_models(z);V<-make_V(ff);rr<-list();gg<-list()
  for(a in seq(40,80,2.5)){
   ref[,age:=a];m<-predict_joint(ff,ref,'Male');f<-predict_joint(ff,ref,'Female');g<-list(est=m$est-f$est,G=m$G-f$G);gg[[as.character(a)]]<-g
   rr[[length(rr)+1]]<-data.table(spec=spec,model=label,analysis='age',age=a,outcome=names(g$est),estimate=g$est,se=sqrt(pmax(0,rowSums((g$G%*%V)*g$G))))
  }
  G<-gg[['75']]$G-gg[['50']]$G;est<-gg[['75']]$est-gg[['50']]$est
  rr[[length(rr)+1]]<-data.table(spec=spec,model=label,analysis='age75_minus50',age=NA_real_,outcome=names(est),estimate=est,se=sqrt(pmax(0,rowSums((G%*%V)*G))))
  fwrite(rbindlist(rr),fn)
  diagnostics<-rbindlist(lapply(names(ff),function(k)data.table(spec=spec,model=label,cause=k,converged=ff[[k]]$converged,max_abs_coefficient=max(abs(ff[[k]]$b)),max_coefficient_se=max(sqrt(diag(ff[[k]]$V))))))
  fwrite(diagnostics,sub('[.]csv$','_checks.csv',fn))
  cat(label,spec,'seconds',proc.time()[3]-tm,'\n');flush.console()
 }
}
