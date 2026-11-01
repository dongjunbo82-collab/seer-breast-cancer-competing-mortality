# Set SEER_ANALYSIS_ROOT to the authorized local SEER input project.
source('work/revision_20261019/engine.R')
source('work/cbc_chain_extension/predict_fast.R')
setDTthreads(1)
w<-'work/jjco_clinical_final'
h<-names(fread(file.path(source_root,'work/cod_extended.csv.gz'),nrows=0));writeLines(h,file.path(w,'available_columns.txt'))
z<-load_data()
fwrite(z[,.(N=.N,deaths=sum(months<=60 & detail %in% causes),breast=sum(months<=60 & detail=='breast'),cardiovascular=sum(months<=60 & detail=='cardiovascular')),by=.(stage,sex)],file.path(w,'stage_event_support.csv'))
rr<-list()
for(c in 1:4)for(i in 1:10){
 label<-sprintf('smc_c%d_i%02d',c,i)
 st<-readRDS(file.path(outroot,'private_imputations',paste0(label,'.rds')))
 zz<-copy(z);stopifnot(identical(zz$row,st$row));for(v in c('stage','grade','er'))zz[,(v):=st[[v]]]
 ref<-zz[sex=='Male',.N,by=vars];ref[,weight:=N/sum(N)]
 ff<-readRDS(file.path(outroot,'fits',paste0(label,'.rds')))$fits
 V<-make_V(ff);m<-predict_joint_fast(ff,ref,'Male');f<-predict_joint_fast(ff,ref,'Female')
 for(ss in c('Male','Female','Difference')){
  a<-switch(ss,Male=m,Female=f,Difference=list(est=m$est-f$est,G=m$G-f$G))
  ks<-c('cardiovascular','other_cancer','respiratory','remaining_other')
  g<-colSums(a$G[ks,,drop=FALSE]);est<-sum(a$est[ks]);stopifnot(abs(est-(a$est['other_total']-a$est['covid']))<1e-12)
  rr[[length(rr)+1]]<-data.table(model=label,analysis='overall_nonCOVID',period='2015-2018',sex=ss,outcome='other_nonCOVID',estimate=est,se=sqrt(max(0,drop(g%*%V%*%g))))
 }
 cat(label,'complete\n');flush.console()
}
x<-rbindlist(rr);fwrite(x,file.path(w,'nonCOVID_per_imputation.csv'))
p<-x[, {B<-var(estimate);T<-mean(se^2)+(1+1/.N)*B;lam<-(1+1/.N)*B/T;df<-if(B<1e-18)Inf else (.N-1)/lam^2;list(m=.N,estimate=mean(estimate),se=sqrt(T),lower=mean(estimate)-qt(.975,df)*sqrt(T),upper=mean(estimate)+qt(.975,df)*sqrt(T))},by=.(analysis,period,sex,outcome)]
fwrite(p,file.path(w,'nonCOVID_pooled.csv'));print(p)
