suppressPackageStartupMessages({library(data.table);library(posterior)})
dest<-Sys.getenv('SEER_EXTENSION_OUT','outputs/chain_stability');args<-commandArgs(TRUE);last<-if(length(args))as.integer(args[1])else 120L
x<-rbindlist(lapply(1:4,function(k)fread(file.path(dest,paste0('estimates_c',k,'.csv')))))[iter<=last]
stopifnot(x[,uniqueN(chain)]==4,all(x[,uniqueN(chain),by=iter]$V1==4));nr<-x[,uniqueN(iter),by=chain];stopifnot(uniqueN(nr$V1)==1);m<-x[,uniqueN(paste(chain,iter))]
pooled<-x[,{
 mat<-as.matrix(dcast(.SD,iter~chain,value.var='estimate')[,-1]);stopifnot(all(is.finite(mat)));B<-var(as.vector(mat));U<-mean(se^2);T<-U+(1+1/m)*B;lambda<-if(T>0)(1+1/m)*B/T else NA_real_;df<-if(B<1e-18)Inf else (m-1)/lambda^2;ess<-posterior::ess_basic(mat,split=FALSE);ac<-apply(mat,2,function(v)as.numeric(acf(v,lag.max=1,plot=FALSE)$acf[2]));mcse<-max(sqrt(B/ess),sqrt(B/m));cm<-colMeans(mat)
 list(m=m,retained_per_chain=nrow(mat),estimate=mean(mat),se=sqrt(T),lower=mean(mat)-qt(.975,df)*sqrt(T),upper=mean(mat)+qt(.975,df)*sqrt(T),within_variance=U,between_variance=B,split_rank_Rhat=posterior::rhat(mat),ESS_mean=ess,ESS_bulk=posterior::ess_bulk(mat),ESS_tail=posterior::ess_tail(mat),lag1_min=min(ac),lag1_max=max(ac),MCSE=mcse,MCSE_SE=mcse/sqrt(T),chain_means_MCSE=sd(cm)/2,chain_mean_min=min(cm),chain_mean_max=max(cm))
},by=.(analysis,age,outcome)]
old<-fread(file.path(Sys.getenv('SEER_REVISION_OUT','outputs/revision_20261019'),'revised_MI.csv'))[sex=='Difference' & period=='2015-2018' & (analysis %in% c('overall','age75_minus50')|analysis=='age'&age%in%c(50,75)),.(analysis,age,outcome,original_estimate=estimate,original_se=se,original_lower=lower,original_upper=upper)]
cmp<-merge(pooled,old,by=c('analysis','age','outcome'));stopifnot(nrow(cmp)==nrow(pooled));cmp[,`:=`(estimate_change_pp=100*(estimate-original_estimate),change_relative_original_SE=(estimate-original_estimate)/original_se,SE_ratio=se/original_se)]
fwrite(cmp,file.path(dest,'extended_chain_comparison.csv'));fwrite(x,file.path(dest,'extended_chain_per_completion.csv'))
cm<-x[,.(estimate_mean=mean(estimate),within_variance=mean(se^2),between_variance=var(estimate)),by=.(chain,analysis,age,outcome)];fwrite(cm,file.path(dest,'extended_chain_means.csv'))
tr<-rbindlist(lapply(1:4,function(k)fread(file.path(dest,paste0('trace_c',k,'.csv')))))[iter>=61 & iter<=last]
cd<-tr[,{
 mat<-as.matrix(dcast(.SD,iter~chain,value.var='proportion')[,-1]);list(n_per_chain=nrow(mat),split_rank_Rhat=posterior::rhat(mat),ESS_bulk=posterior::ess_bulk(mat))
},by=.(variable,sex,category)];fwrite(cd,file.path(dest,'extended_category_diagnostics.csv'));writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
print(cmp[analysis=='overall'&outcome%in%c('all','other_total','cardiovascular')|analysis=='age75_minus50'&outcome%in%c('other_total','cardiovascular')|analysis=='age'&age==50&outcome=='cardiovascular',.(analysis,age,outcome,m,estimate,se,split_rank_Rhat,ESS_bulk,MCSE_SE,estimate_change_pp,SE_ratio)])
cat('MAX SCALAR RHAT',max(cmp$split_rank_Rhat),'MAX CATEGORY RHAT',max(cd$split_rank_Rhat),'\n')
