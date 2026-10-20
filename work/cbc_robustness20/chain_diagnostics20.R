suppressPackageStartupMessages({library(data.table);library(posterior)})
dest<-'outputs/cbc_robustness20';dir.create(dest,recursive=TRUE,showWarnings=FALSE)
input<-file.path(Sys.getenv('SEER_REVISION_OUT','outputs/revision_20261019'),'revised_MI_individual.csv')
if(!file.exists(input))input<-'outputs/tables/revision19_revised_MI_individual.csv'
raw<-fread(input)

raw[,chain:=as.integer(sub('^smc_c([1-4]).*','\\1',model))]
raw[,retained:=as.integer(sub('^smc_c[1-4]_i([0-9]+).*','\\1',model))]
# Ten retained draws per chain map to sweeps 22,24,...,40.
x<-raw[period=='2015-2018' & sex=='Difference' & analysis %in% c('overall','age75_minus50')]
x<-rbind(x,raw[analysis=='age' & sex=='Difference' & age %in% c(50,75) & outcome %in% c('other_total','cardiovascular')],raw[analysis=='overall' & sex=='Female' & period=='2010-2014' & outcome=='breast'])
ans<-x[,{
 mat<-as.matrix(dcast(.SD,retained~chain,value.var='estimate')[,-1]);stopifnot(all(dim(mat)==c(10,4)),all(is.finite(mat)))
 ac<-apply(mat,2,function(v)as.numeric(acf(v,lag.max=1,plot=FALSE)$acf[2]))
 ess_split<-posterior::ess_mean(mat);ess<-posterior::ess_basic(mat,split=FALSE)
 mcse<-sqrt(var(as.vector(mat))/ess)
 B<-var(as.vector(mat));U<-mean(se^2);T<-U+(1+1/40)*B;iid<-sqrt(B/40)
 # With only ten retained draws/chain, do not claim a precision gain from estimated negative correlation.
 safe<-max(mcse,iid)
 list(n=40,chains=4,retained_per_chain=10,sweeps='22,24,...,40',lag1_min=min(ac),lag1_max=max(ac),lag1_mean=mean(ac),ESS_unsplit=ess,ESS_split=ess_split,split_rank_Rhat=posterior::rhat(mat),iid_MCSE=iid,autocorrelation_MCSE=mcse,conservative_MCSE=safe,total_SE=sqrt(T),autocorrelation_MCSE_SE=mcse/sqrt(T),conservative_MCSE_SE=safe/sqrt(T))
},by=.(analysis,age,sex,period,outcome)]
fwrite(ans,file.path(dest,'chain_diagnostics_key_estimands.csv'))
# Independent chains: estimate MCSE from variability of the four chain means as a coarse cross-check.
cm<-x[,.(mean=mean(estimate)),by=.(analysis,age,sex,period,outcome,chain)]
check<-cm[,.(chain_means_MCSE=sd(mean)/sqrt(.N)),by=.(analysis,age,sex,period,outcome)]
fwrite(check,file.path(dest,'chain_means_mcse_check.csv'))
writeLines(capture.output(sessionInfo()),file.path(dest,'diagnostic_sessionInfo.txt'))
print(ans[analysis=='age75_minus50' & outcome %in% c('other_total','cardiovascular') | analysis=='overall' & sex=='Difference' & outcome %in% c('all','other_total')])
