suppressPackageStartupMessages(library(data.table))
dest<-'outputs/cbc_robustness20'
fs<-list.files(dest,pattern='^(two_knots|shifted_knots)_smc_c[1-4]_i[0-9]+[.]csv$',full.names=TRUE);stopifnot(length(fs)==80)
x<-rbindlist(lapply(fs,fread));stopifnot(all(x[, .N,by=.(spec,analysis,age,outcome)]$N==40))
base<-fread(file.path(Sys.getenv('SEER_REVISION_OUT','outputs/revision_20261019'),'revised_MI_individual.csv'))[period=='2015-2018' & sex=='Difference' & analysis %in% c('age','age75_minus50'),.(spec='primary_knots',model,analysis,age,outcome,estimate,se)]
x<-rbind(x,base);fwrite(x,file.path(dest,'spline_sensitivity_per_imputation.csv'))
ans<-x[,{
 n<-.N;Q<-mean(estimate);U<-mean(se^2);B<-var(estimate);T<-U+(1+1/n)*B;lam<-(1+1/n)*B/T;df<-if(B<1e-18)Inf else (n-1)/lam^2
 list(m=n,estimate=Q,se=sqrt(T),lower=Q-qt(.975,df)*sqrt(T),upper=Q+qt(.975,df)*sqrt(T))
},by=.(spec,analysis,age,outcome)]
fwrite(ans,file.path(dest,'spline_sensitivity_pooled.csv'))
ck<-rbindlist(lapply(list.files(dest,pattern='_checks.csv$',full.names=TRUE),fread));stopifnot(nrow(ck)==560,all(ck$converged),all(is.finite(ck$max_abs_coefficient)),all(is.finite(ck$max_coefficient_se)))
fwrite(ck,file.path(dest,'spline_fit_diagnostics.csv'))
fwrite(ck[,.(n_fits=.N,all_converged=all(converged),max_abs_coefficient=max(max_abs_coefficient),max_coefficient_se=max(max_coefficient_se)),by=.(spec,cause)],file.path(dest,'spline_fit_summary.csv'))
print(ans[analysis=='age75_minus50' & outcome %in% c('other_total','cardiovascular')])
writeLines(capture.output(sessionInfo()),file.path(dest,'spline_sessionInfo.txt'))
