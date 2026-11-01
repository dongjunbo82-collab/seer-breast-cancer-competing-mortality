# Set SEER_ANALYSIS_ROOT to your authorized local input project.

source('work/revision_20261019/engine.R');library(survival);setDTthreads(1)
w<-'work/jjco_recent_cv';dir.create(w,showWarnings=FALSE)
d<-load_data()[year>=2015];d[,age10:=(age-65)/10];d[,age_group:=cut(age,c(17,54,64,74,Inf),labels=c('18-54','55-64','65-74','75+'))]
fwrite(d[,.(N=.N,cardiovascular=sum(detail=='cardiovascular'&months<=60),all_deaths=sum(detail%in%causes&months<=60)),by=.(sex,age_group)],file.path(w,'support.csv'))
ks<-c('breast','cardiovascular','remaining_deaths');d[,endpoint:=fifelse(detail=='breast','breast',fifelse(detail=='cardiovascular','cardiovascular',fifelse(detail%in%causes,'remaining_deaths','censor')))];d[time>5,endpoint:='censor'];d[,time:=pmin(time,5)]
# Observed Aalen-Johansen estimates, retaining all other deaths as competing events.
aj<-list()
for(ss in levels(d$sex))for(gg in levels(d$age_group)){
 z<-d[sex==ss&age_group==gg];fit<-survfit(Surv(time,factor(endpoint,levels=c('censor',ks)))~1,data=z)
 j<-length(fit$time);k<-match('cardiovascular',fit$states);est<-fit$pstate[j,k];se<-fit$std.err[j,k]
 if(est>0&&est<1){q<-log(-log1p(-est));s<-se/((1-est)*(-log1p(-est)));lo<--expm1(-exp(q-1.96*s));hi<--expm1(-exp(q+1.96*s))}else{lo<-hi<-NA_real_}
 aj[[length(aj)+1]]<-data.table(sex=ss,age_group=gg,N=nrow(z),cardiovascular_events=sum(z$endpoint=='cardiovascular'),estimate=est,se=se,lower=lo,upper=hi)
}
fwrite(rbindlist(aj),file.path(w,'observed_AJ.csv'))
form<-~sex*age10+interval+sex:interval+stage+grade+er+race
predict3<-function(ff,ref,ss,a){
 p<-length(ff[[1]]$b);K<-3*p;n<-nrow(ref);z<-copy(ref);z[,sex:=factor(ss,levels=levels(ref$sex))];z[,age10:=(a-65)/10]
 S<-rep(1,n);F<-matrix(0,n,3);DS<-matrix(0,n,K);DF<-lapply(1:3,function(k)matrix(0,n,K))
 for(j in 1:5){z[,interval:=factor(j,levels=1:5)];X<-model.matrix(form,z);h<-sapply(ff,function(f)exp(drop(X%*%f$b)));H<-rowSums(h);q<-exp(-H);v<--expm1(-H)/H;vp<-(q*H+expm1(-H))/H^2;DH<-do.call(cbind,lapply(1:3,function(k)X*h[,k]))
 for(k in 1:3){DF[[k]]<-DF[[k]]+DS*(v*h[,k])+DH*(S*vp*h[,k]);ix<-((k-1)*p+1):(k*p);DF[[k]][,ix]<-DF[[k]][,ix]+X*(S*v*h[,k]);F[,k]<-F[,k]+S*v*h[,k]};DS<-DS*q-DH*(S*q);S<-S*q}
 list(est=setNames(colSums(F*ref$weight),ks),G=t(sapply(DF,function(g)colSums(g*ref$weight))))
}
rr<-list();diag<-list()
for(mode in c('unknown_category','complete_case')){
 z<-copy(d);if(mode=='complete_case')z<-z[stage!='Unknown'&grade!='Unknown'&er!='Unknown'];z<-droplevels(z)
 pt<-rbindlist(lapply(1:5,function(j){a<-z[time>j-1];a[,exposure:=pmin(time,j)-(j-1)];for(k in ks)a[,(k):=as.integer(endpoint==k&time<=j)];a<-a[,c(list(exposure=sum(exposure)),lapply(.SD,sum)),by=.(sex,age10,stage,grade,er,race),.SDcols=ks];a[,interval:=factor(j,levels=1:5)];a}))
 X<-model.matrix(form,pt);ff<-list();p<-ncol(X);V<-matrix(0,3*p,3*p)
 for(k in 1:3){fit<-glm.fit(X,pt[[ks[k]]],offset=log(pt$exposure),family=poisson(),control=glm.control(maxit=100));class(fit)<-c('glm','lm');vv<-vcov(fit);stopifnot(fit$converged,all(is.finite(coef(fit))),all(is.finite(vv)));ff[[k]]<-list(b=coef(fit),V=vv);ix<-((k-1)*p+1):(k*p);V[ix,ix]<-vv;diag[[length(diag)+1]]<-data.table(model=mode,cause=ks[k],converged=fit$converged,events=sum(pt[[ks[k]]]),n_male=sum(z$sex=='Male'),n_female=sum(z$sex=='Female'),max_coefficient_se=max(sqrt(diag(vv))))}
 ref<-z[sex=='Male',.N,by=.(sex,stage,grade,er,race)];ref[,weight:=N/sum(N)];ages<-list()
 for(a in seq(40,80,2.5)){m<-predict3(ff,ref,'Male',a);f<-predict3(ff,ref,'Female',a);g<-list(est=m$est-f$est,G=m$G-f$G);ages[[as.character(a)]]<-g
 for(ss in c('Male','Female','Difference')){o<-switch(ss,Male=m,Female=f,Difference=g);se<-sqrt(pmax(0,rowSums((o$G%*%V)*o$G)));rr[[length(rr)+1]]<-data.table(model=mode,analysis='age',age=a,sex=ss,outcome=ks,estimate=o$est,se=se,lower=o$est-1.96*se,upper=o$est+1.96*se)}}
 g<-list(est=ages[['75']]$est-ages[['50']]$est,G=ages[['75']]$G-ages[['50']]$G);se<-sqrt(pmax(0,rowSums((g$G%*%V)*g$G)));rr[[length(rr)+1]]<-data.table(model=mode,analysis='age75_minus50',age=NA_real_,sex='Difference',outcome=ks,estimate=g$est,se=se,lower=g$est-1.96*se,upper=g$est+1.96*se)
 # Numerical gradient check at age65 for a cardiovascular coefficient.
 o<-predict3(ff,ref,'Male',65);j<-which(names(ff[[2]]$b)=='sexMale:age10');fp<-fm<-ff;fp[[2]]$b[j]<-fp[[2]]$b[j]+1e-5;fm[[2]]$b[j]<-fm[[2]]$b[j]-1e-5;fd<-(predict3(fp,ref,'Male',65)$est-predict3(fm,ref,'Male',65)$est)/2e-5;stopifnot(max(abs(fd-o$G[,p+j]))<1e-7)
 saveRDS(list(fits=ff,reference=ref),file.path(w,paste0(mode,'_private_fit.rds')));cat(mode,'completed\n')
}
fwrite(rbindlist(rr),file.path(w,'model_estimates.csv'));fwrite(rbindlist(diag),file.path(w,'fit_diagnostics.csv'));print(rbindlist(rr)[outcome=='cardiovascular'&(analysis=='age75_minus50'|age%in%c(50,65,75))&sex=='Difference'])
