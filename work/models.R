suppressPackageStartupMessages({library(data.table);library(splines);library(survival)})
setDTthreads(4); set.seed(20261006)
dir.create('work/models',showWarnings=FALSE)
d<-fread('work/extended.csv.gz')
d[,period:=cut(year,c(1999,2004,2009,2014,2018,2020,2023),labels=c('2000-2004','2005-2009','2010-2014','2015-2018','2019-2020','2021-2023'))]
d[,time:=pmax(months,0.5)/12]
for(v in c('sex','period','race','stage','grade','er'))d[,(v):=factor(get(v))]
d[,sex:=relevel(sex,'Female')];d[,stage:=relevel(stage,'Localized')]
fwrite(d[,.N,by=.(sex,race)],'work/race_check.csv')
main<-d[mainhist==1 & year<=2018];main[,period:=droplevels(period)]
saveRDS(main,'work/main.rds',compress=FALSE)
vars<-c('sex','period','age','race','stage','grade','er')
make_person_time<-function(z,horizon=5){
 rbindlist(lapply(seq_len(ceiling(horizon)),function(j){
  a<-z[time>(j-1)];a[,exposure:=pmin(time,j,horizon)-(j-1)]
  a[,`:=`(bc=as.integer(event==1 & time<=min(j,horizon)),other=as.integer(event==2 & time<=min(j,horizon)),unknown=as.integer(event==3 & time<=min(j,horizon)))]
  a<-a[,.(exposure=sum(exposure),bc=sum(bc),other=sum(other),unknown=sum(unknown)),by=vars]
  a[,interval:=factor(j,levels=seq_len(ceiling(horizon)))];a
 }))
}
base<- '~ sex*period + ns(age,knots=c(50,65,75),Boundary.knots=c(18,85)) + interval + sex:interval'
forms<-list(age=as.formula(base),age_stage=as.formula(paste(base,'+ stage + stage:interval')),full=as.formula(paste(base,'+ stage + stage:interval + grade + er + race')))
fit_hazards<-function(z,form,horizon=5){
 pt<-make_person_time(z,horizon);cat('person-time rows',nrow(pt),'\n');flush.console()
 X<-model.matrix(form,pt); off<-log(pt$exposure); ans<-list()
 for(outcome in c('bc','other','unknown')){
  f<-glm.fit(X,pt[[outcome]],offset=off,family=poisson(),control=glm.control(maxit=40))
  class(f)<-c('glm','lm'); V<-vcov(f)
  stopifnot(f$converged,!anyNA(coef(f)),all(is.finite(V)))
  ans[[outcome]]<-list(b=coef(f),V=V,formula=form,terms=terms(form),contrasts=attr(X,'contrasts'),deviance=f$deviance,df=f$df.residual,converged=f$converged)
  cat('fit',outcome,'deviance/df',f$deviance/f$df.residual,'\n');flush.console()
 }
 ans
}
predict_risk<-function(fits,ref,sex0,period0,horizon=5){
 n<-nrow(ref);p<-length(fits[[1]]$b);K<-3*p
 zz<-copy(ref);zz[,sex:=factor(sex0,levels=levels(ref$sex))];zz[,period:=factor(period0,levels=levels(ref$period))]
 S<-rep(1,n);F<-matrix(0,n,3);DS<-matrix(0,n,K);DF<-lapply(1:3,function(k)matrix(0,n,K))
 for(j in seq_len(ceiling(horizon))){
  zz[,interval:=factor(j,levels=seq_len(ceiling(horizon)))];X<-model.matrix(fits[[1]]$formula,zz,contrasts.arg=fits[[1]]$contrasts)
  hh<-sapply(fits,function(f)exp(drop(X%*%f$b)))*min(1,horizon-j+1);H<-rowSums(hh);q<-exp(-H);a<--expm1(-H)/H;ap<-(q*H+expm1(-H))/H^2
  DH<-do.call(cbind,lapply(1:3,function(k)X*hh[,k]))
  for(k in 1:3){DF[[k]]<-DF[[k]]+DS*(a*hh[,k])+DH*(S*ap*hh[,k]);ii<-((k-1)*p+1):(k*p);DF[[k]][,ii]<-DF[[k]][,ii]+X*(S*a*hh[,k]);F[,k]<-F[,k]+S*a*hh[,k]}
  DS<-DS*q-DH*(S*q);S<-S*q
 }
 w<-ref$weight;est<-c(all=sum((1-S)*w),bc=sum(F[,1]*w),other=sum(F[,2]*w),unknown=sum(F[,3]*w))
 G<-rbind(-colSums(DS*w),t(sapply(DF,function(g)colSums(g*w))))
 stopifnot(abs(est[1]-sum(est[-1]))<1e-9)
 list(est=est,G=G)
}
analyse<-function(z,label,form=forms$full,horizon=5,reference=NULL){
 cat('START',label,nrow(z),'\n');flush.console()
 z<-droplevels(copy(z));z[,period:=droplevels(period)]
 fitfile<-paste0('work/models/',label,'_fit.rds')
 ff<-if(file.exists(fitfile))readRDS(fitfile) else {v<-fit_hazards(z,form,horizon);saveRDS(v,fitfile);v}
 ref<-if(is.null(reference)) z[sex=='Male',.N,by=vars] else copy(reference)
 ref[,weight:=N/sum(N)];ref[,period:=factor(as.character(period),levels=levels(z$period))]
 p<-length(ff[[1]]$b);V<-matrix(0,3*p,3*p)
 for(k in 1:3){ii<-((k-1)*p+1):(k*p);V[ii,ii]<-ff[[k]]$V}
 tab<-list();prs<-list();gps<-list();per<-levels(z$period)
 getrows<-function(obj,kind,period,sex='Difference'){
  se<-sqrt(pmax(0,rowSums((obj$G%*%V)*obj$G)))
  data.table(model=label,kind=kind,period=period,sex=sex,outcome=names(obj$est),estimate=obj$est,se=se,lower=obj$est-1.96*se,upper=obj$est+1.96*se,p=2*pnorm(-abs(obj$est/se)))
 }
 for(pp in per){
  for(ss in c('Female','Male')){o<-predict_risk(ff,ref,ss,pp,horizon);prs[[paste(ss,pp)]]<-o;tab[[length(tab)+1]]<-getrows(o,'risk',pp,ss)}
  m<-prs[[paste('Male',pp)]];f<-prs[[paste('Female',pp)]];g<-list(est=m$est-f$est,G=m$G-f$G);gps[[pp]]<-g
  tab[[length(tab)+1]]<-getrows(g,'gap',pp)
 }
 first<-gps[[1]];last<-gps[[length(gps)]];chg<-list(est=last$est-first$est,G=last$G-first$G)
 tab[[length(tab)+1]]<-getrows(chg,'change',paste(tail(per,1),'-',per[1]))
 tt<-rbindlist(tab);fwrite(tt,paste0('work/models/',label,'_results.csv'))
 # Three-degree-of-freedom test of whether standardized sex gaps vary by period.
 joint<-rbindlist(lapply(1:4,function(o){J<-do.call(rbind,lapply(gps[-1],function(g)g$G[o,]-first$G[o,]));b<-sapply(gps[-1],function(g)g$est[o]-first$est[o]);W<-drop(t(b)%*%solve(J%*%V%*%t(J),b));data.table(model=label,outcome=names(first$est)[o],chisq=W,df=length(b),p=pchisq(W,length(b),lower.tail=FALSE))}))
 fwrite(joint,paste0('work/models/',label,'_joint.csv'))
 saveRDS(list(ref=ref,predictions=prs,gaps=gps,V=V),paste0('work/models/',label,'_prediction.rds'))
 print(tt[kind=='change']);cat('DONE',label,'\n');flush.console();invisible(tt)
}
if(!identical(Sys.getenv('LOAD_ONLY'),'1')){
 for(mm in names(forms))analyse(main,mm,forms[[mm]])
}
