# Secondary national background benchmark; not a relative-survival or attributable-risk analysis.
# Set SEER_ANALYSIS_ROOT to the locally authorized input project before running.
source('work/revision_20261019/engine.R');source('work/jjco_value_revision/predict_fast.R')
setDTthreads(1);dest<-'work/jjco_value_revision';dir.create(dest,showWarnings=FALSE)
z<-load_data();N<-z[sex=='Male',.N];eligible<-z[sex=='Male'&age<85,.N];excluded<-N-eligible
agew<-z[sex=='Male'&age<85,.N,by=age][,weight:=N/sum(N)]
yw<-z[sex=='Male'&age<85&year>=2015,.N,by=year][,weight:=N/sum(N)]
q<-fread('outputs/JJCO_Submission/Supplementary_Data_and_Code/outputs/tables/life_table_qx.csv');setkey(q,sex,year,age)
a<-CJ(age=agew$age,year=yw$year);a<-merge(a,agew[,.(age,age_weight=weight)],by='age');a<-merge(a,yw[,.(year,year_weight=weight)],by='year');a[,weight:=age_weight*year_weight]
bg<-rbindlist(lapply(c(0,.5),function(offset)rbindlist(lapply(c('Male','Female'),function(ss){
 S<-rep(1,nrow(a));for(j in 0:4){q0<-q[data.table(sex=ss,year=a$year+j,age=a$age+j)]$q;q1<-q[data.table(sex=ss,year=a$year+j+1,age=a$age+j)]$q;stopifnot(!anyNA(q0),!anyNA(q1));S<-S*(1-q0)^(1-offset)*(1-q1)^offset}
 data.table(offset=offset,sex=ss,estimate=sum(a$weight*(1-S)))
}))))
fwrite(bg,file.path(dest,'background_expected.csv'));fwrite(agew,file.path(dest,'background_age_weights.csv'));fwrite(yw,file.path(dest,'background_year_weights.csv'))
rr<-list();check<-NULL
for(c in 1:4)for(i in 1:10){
 label<-sprintf('smc_c%d_i%02d',c,i);st<-readRDS(file.path(outroot,'private_imputations',paste0(label,'.rds')));stopifnot(identical(st$row,z$row));for(v in c('stage','grade','er'))z[,(v):=st[[v]]]
 ff<-readRDS(file.path(outroot,'fits',paste0(label,'.rds')))$fits
 ref<-z[sex=='Male'&age<85,.N,by=vars];ref[,weight:=N/sum(N)]
 if(c==1&&i==1){
  test<-copy(ref[seq_len(min(100,.N))]);test[,weight:=weight/sum(weight)];slow<-predict_joint(ff,test,'Male');fast<-predict_joint_fast(ff,test,'Male');check<-data.table(check=c('estimate','gradient'),max_error=c(max(abs(slow$est-fast$est)),max(abs(slow$G-fast$G))));stopifnot(max(check$max_error)<1e-9);fwrite(check,file.path(dest,'background_gradient_validation.csv'))
 }
 m<-predict_joint_fast(ff,ref,'Male');f<-predict_joint_fast(ff,ref,'Female');V<-make_V(ff)
 for(ss in c('Male','Female','Difference')){ans<-if(ss=='Male')m else if(ss=='Female')f else list(est=m$est-f$est,G=m$G-f$G);rr[[length(rr)+1]]<-data.table(model=label,sex=ss,outcome=names(ans$est),estimate=ans$est,se=sqrt(pmax(0,rowSums((ans$G%*%V)*ans$G))))}
 fwrite(rbindlist(rr),file.path(dest,'background_restricted_per_imputation.csv'));cat(label,'reference profiles',nrow(ref),'\n');flush.console()
}
x<-rbindlist(rr);pooled<-x[,{
 B<-var(estimate);U<-mean(se^2);T<-U+(1+1/.N)*B;lam<-if(T>0)(1+1/.N)*B/T else 0;df<-if(B<1e-18)Inf else (.N-1)/lam^2;list(m=.N,estimate=mean(estimate),se=sqrt(T),lower=mean(estimate)-qt(.975,df)*sqrt(T),upper=mean(estimate)+qt(.975,df)*sqrt(T))
},by=.(sex,outcome)]
fwrite(pooled,file.path(dest,'background_restricted_pooled.csv'));fwrite(data.table(male_total=N,male_age18_84=eligible,male_85plus_excluded=excluded),file.path(dest,'background_eligibility.csv'))
print(pooled[outcome%in%c('all','breast','other_total')]);print(bg);writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
