source('work/revision_20261019/engine.R')
files<-list.files(file.path(outroot,'fits'),pattern='^smc_c[1-4]_i[0-9]+[.]rds$',full.names=TRUE)
stopifnot(length(files)==40);objs<-lapply(files,readRDS)
pool<-function(x,groups)x[,{
 m<-.N;Q<-mean(estimate);U<-mean(se^2);B<-var(estimate);T<-U+(1+1/m)*B;lam<-if(T>0)(1+1/m)*B/T else 0
 df<-if(B<1e-18)Inf else (m-1)/lam^2
 list(estimate=Q,se=sqrt(T),lower=Q-qt(.975,df)*sqrt(T),upper=Q+qt(.975,df)*sqrt(T),p=if(T>0)2*pt(-abs(Q/sqrt(T)),df) else NA_real_,m=m,FMI=lam+(1-lam)*2/(df+3),MCSE_SE=if(T>0)sqrt(B/m/T) else 0)
},by=groups]
raw<-rbindlist(lapply(objs,`[[`,'results'))
fwrite(raw,file.path(outroot,'revised_MI_individual.csv'))
ans<-pool(raw,c('analysis','age','sex','period','outcome'))
ans[,p_holm:=NA_real_]
for(a in c('age75_minus50','supported_age75_minus50','latest_supported_age75_minus50')){
 for(oc in list(c('all','breast','other_total','unknown'),causes[2:6])){
  ii<-which(ans$analysis==a & ans$outcome %in% oc);ans[ii,p_holm:=p.adjust(p,'holm')]
 }
}
ans[analysis=='cause_comparison',p_holm:=p.adjust(p,'holm'),by=age]
fwrite(ans,file.path(outroot,'revised_MI.csv'))
ss<-rbindlist(lapply(c('unknown_category','complete_case'),function(lab){
 o<-readRDS(file.path(outroot,'fits',paste0(lab,'.rds')));x<-copy(o$results)
 x[,`:=`(lower=estimate-1.96*se,upper=estimate+1.96*se,p=2*pnorm(-abs(estimate/se)))];x[,p_holm:=NA_real_]
 for(a in c('age75_minus50','supported_age75_minus50'))for(oc in list(c('all','breast','other_total','unknown'),causes[2:6])){
  ii<-which(x$analysis==a & x$outcome %in% oc);x[ii,p_holm:=p.adjust(p,'holm')]
 };x
}))
fwrite(ss,file.path(outroot,'revised_sensitivity.csv'))
su<-rbindlist(lapply(objs,`[[`,'support'));fwrite(su,file.path(outroot,'support_individual.csv'))
fwrite(su[,lapply(.SD,mean),by=age,.SDcols=c('reference_weight_any_missing','reference_weight_below5','latest_zero_weight','latest_below5_weight','male_events')],file.path(outroot,'support_summary.csv'))
meta<-rbindlist(lapply(seq_along(objs),function(i){o<-objs[[i]];data.table(imputation=basename(files[i]),n=o$n,supported_reference_weight=o$supported_reference_weight,latest_supported_reference_weight=o$latest_supported_reference_weight,gradient_error=o$gradient_error,converged=all(sapply(o$fits,`[[`,'converged')))}))
fwrite(meta,file.path(outroot,'revised_checks.csv'))
ix<-rbindlist(lapply(seq_along(objs),function(i)copy(objs[[i]]$intersection_counts)[,imputation:=i]));fwrite(ix,file.path(outroot,'intersection_counts.csv'))
tr<-rbindlist(lapply(1:4,function(k)fread(file.path(outroot,paste0('trace_c',k,'.csv')))))
# Split-chain Rhat on missing-category proportions after the first ten sweeps.
# correlation of retained summaries. These diagnose chains, not MAR validity.
dg<-tr[iter>=11,{
 mat<-dcast(.SD,iter~chain,value.var='proportion');M<-as.matrix(mat[,-1]);half<-floor(nrow(M)/2);M<-cbind(M[seq_len(half),,drop=FALSE],M[half+seq_len(half),,drop=FALSE]);n<-nrow(M)
 W<-mean(apply(M,2,var));B<-n*var(colMeans(M));rh<-if(W>0)sqrt(((n-1)/n*W+B/n)/W) else 1
 list(split_Rhat=rh)
},by=.(variable,sex,category)]
fwrite(dg,file.path(outroot,'imputation_diagnostics.csv'))
cf<-list.files(file.path(outroot,'fits'),pattern='_model_comparison.csv$',full.names=TRUE)
if(length(cf)==40){
 cc<-rbindlist(lapply(cf,fread));fwrite(pool(cc,c('model','analysis','sex','age','outcome')),file.path(outroot,'model_comparison.csv'))
}
stopifnot(all(meta$converged),max(meta$gradient_error)<1e-6,all(ans$m==40))
co<-ans[analysis %in% c('overall','age'),.(err=estimate[outcome=='all']-sum(estimate[outcome %in% causes])),by=.(analysis,age,sex,period)]
stopifnot(max(abs(co$err))<1e-10)
print(ans[analysis=='overall' & period=='2015-2018' & sex=='Difference'])
print(ans[analysis=='age75_minus50']);print(dg);print(meta[,.(min_support=min(supported_reference_weight),max_grad=max(gradient_error))])
writeLines(capture.output(sessionInfo()),file.path(outroot,'sessionInfo.txt'))
