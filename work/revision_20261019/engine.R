suppressPackageStartupMessages({library(data.table);library(splines);library(nnet);library(MASS)})
setDTthreads(2)
source_root<-Sys.getenv('SEER_ANALYSIS_ROOT','.')
outroot<-Sys.getenv('SEER_REVISION_OUT','outputs/revision_20261019')
dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
causes<-c('breast','cardiovascular','other_cancer','respiratory','covid','remaining_other','unknown')
vars<-c('sex','period','age','race','stage','grade','er')
form<-~sex*period+sex*ns(age,knots=c(50,65,75),Boundary.knots=c(18,85))+I((age-65)*(sex=='Male'))+interval+sex:interval+stage+stage:interval+grade+er+race
load_data<-function(){
 d<-fread(file.path(source_root,'work/cod_extended.csv.gz'))[mainhist==1 & year<=2018]
 old<-readRDS(file.path(source_root,'work/main.rds'));stopifnot(identical(d$row,old$row));rm(old)
 d[,period:=factor(cut(year,c(1999,2004,2009,2014,2018),labels=c('2000-2004','2005-2009','2010-2014','2015-2018')))]
 d[,time:=pmax(months,.5)/12]
 for(v in c('sex','race','stage','grade','er'))d[,(v):=factor(get(v))]
 d[,sex:=relevel(sex,'Female')];d[,stage:=relevel(stage,'Localized')];d
}
fit_models<-function(d,draw=FALSE){
 pt<-rbindlist(lapply(1:5,function(j){
  a<-d[time>j-1];a[,exposure:=pmin(time,j)-(j-1)]
  for(k in causes)a[,(k):=as.integer(detail==k & time<=j)]
  a<-a[,c(list(exposure=sum(exposure)),lapply(.SD,sum)),by=vars,.SDcols=causes]
  a[,interval:=factor(j,levels=1:5)];a
 }))
 X<-model.matrix(form,pt);p<-ncol(X);ff<-list()
 for(k in causes){
  use<-if(k=='covid')which(pt$period=='2015-2018' & as.integer(pt$interval)>=2) else seq_len(nrow(pt))
  xx<-X[use,,drop=FALSE]
  # Only 11 male COVID deaths and 24 unknown-cause deaths in the full cohort:
  # do not estimate four additional sex-age coefficients for those sparse causes.
  nonlinear<-grepl('sexMale:ns|ns.*:sexMale',colnames(X));linear<-grepl('^I\\(',colnames(X))
  allowed<-if(k %in% c('covid','unknown'))which(!nonlinear & !linear) else if(k=='respiratory')which(!nonlinear) else which(!linear)
  qq<-qr(xx[,allowed,drop=FALSE]);keep<-sort(allowed[qq$pivot[seq_len(qq$rank)]])
  f<-glm.fit(xx[,keep,drop=FALSE],pt[[k]][use],offset=log(pt$exposure[use]),family=poisson(),control=glm.control(maxit=80));class(f)<-c('glm','lm')
  stopifnot(f$converged,all(is.finite(coef(f))),all(is.finite(vcov(f))))
  b<-numeric(p);names(b)<-colnames(X);b[keep]<-if(draw)drop(MASS::mvrnorm(1,coef(f),vcov(f))) else coef(f)
  V<-matrix(0,p,p);V[keep,keep]<-vcov(f)
  ff[[k]]<-list(b=b,V=V,events=sum(pt[[k]]),converged=f$converged,rank=length(keep))
 }
 ff
}
# Exact observed competing-risk likelihood for a fixed completed covariate pattern.
# Enumerating the finite categories avoids rejection-sampling limits.
loglik<-function(ff,z){
 n<-nrow(z);if(!n)return(numeric());zz<-copy(z);ll<-numeric(n)
 for(j in 1:5){
  zz[,interval:=factor(j,levels=1:5)];X<-model.matrix(form,zz)
  eta<-sapply(ff,function(f)drop(X%*%f$b));hh<-exp(eta)
  hh[zz$period!='2015-2018' | j==1,5]<-0
  expo<-pmax(0,pmin(zz$time,j)-(j-1));ll<-ll-rowSums(hh)*expo
  ev<-which(zz$time>j-1 & zz$time<=j & zz$detail %in% causes)
  if(length(ev)){kk<-match(zz$detail[ev],causes);stopifnot(all(hh[cbind(ev,kk)]>0));ll[ev]<-ll[ev]+eta[cbind(ev,kk)]}
 }
 ll
}
# Exact category likelihoods with one shared design calculation per interval.
# The imputed covariates enter as main effects or stage-by-interval effects, so
# changing a category shifts log hazards by a vector constant across patients.
category_loglik<-function(ff,z,variable,categories){
 zz<-copy(z);zz[,(variable):=factor(categories[1],levels=levels(z[[variable]]))]
 n<-nrow(zz);L<-length(categories);ll<-matrix(0,n,L);B<-do.call(cbind,lapply(ff,`[[`,'b'))
 for(j in 1:5){
  zz[,interval:=factor(j,levels=1:5)];X<-model.matrix(form,zz);eta<-X%*%B
  template<-zz[rep(1,L)];template[,(variable):=factor(categories,levels=levels(z[[variable]]))]
  shift<-(model.matrix(form,template)-matrix(X[1,],L,ncol(X),byrow=TRUE))%*%B
  hazard<-exp(eta);hazard[zz$period!='2015-2018' | j==1,5]<-0
  exposure<-pmax(0,pmin(zz$time,j)-(j-1))
  ll<-ll-(hazard*exposure)%*%t(exp(shift))
  ev<-which(zz$time>j-1 & zz$time<=j & zz$detail %in% causes)
  if(length(ev)){
   k<-match(zz$detail[ev],causes)
   ll[ev,]<-ll[ev,,drop=FALSE]+matrix(eta[cbind(ev,k)],length(ev),L)+t(shift[,k,drop=FALSE])
  }
 }
 ll
}
predict_joint<-function(ff,ref,ss,pp='2015-2018'){
 p<-length(ff[[1]]$b);K<-7*p;n<-nrow(ref);z<-copy(ref)
 z[,sex:=factor(ss,levels=levels(ref$sex))];z[,period:=factor(pp,levels=levels(ref$period))]
 S<-rep(1,n);F<-matrix(0,n,7);DS<-matrix(0,n,K);DF<-lapply(1:7,function(k)matrix(0,n,K))
 for(j in 1:5){
  z[,interval:=factor(j,levels=1:5)];X<-model.matrix(form,z)
  hh<-sapply(ff,function(f)exp(drop(X%*%f$b)));if(pp!='2015-2018' || j==1)hh[,5]<-0
  H<-rowSums(hh);q<-exp(-H);a<--expm1(-H)/H;ap<-(q*H+expm1(-H))/H^2
  DH<-do.call(cbind,lapply(1:7,function(k)X*hh[,k]))
  for(k in 1:7){DF[[k]]<-DF[[k]]+DS*(a*hh[,k])+DH*(S*ap*hh[,k]);ii<-((k-1)*p+1):(k*p);DF[[k]][,ii]<-DF[[k]][,ii]+X*(S*a*hh[,k]);F[,k]<-F[,k]+S*a*hh[,k]}
  DS<-DS*q-DH*(S*q);S<-S*q
 }
 w<-ref$weight;est<-c(all=sum((1-S)*w),setNames(colSums(F*w),causes))
 G<-rbind(-colSums(DS*w),t(sapply(DF,function(g)colSums(g*w))));rownames(G)<-names(est)
 est<-c(est,other_total=sum(est[causes[2:6]]));G<-rbind(G,other_total=colSums(G[causes[2:6],,drop=FALSE]))
 stopifnot(abs(est[1]-sum(est[causes]))<1e-10,all(est>=0),all(est<=1));list(est=est,G=G)
}
make_V<-function(ff){p<-length(ff[[1]]$b);V<-matrix(0,7*p,7*p);for(k in 1:7){ix<-((k-1)*p+1):(k*p);V[ix,ix]<-ff[[k]]$V};V}
results<-function(ff,z,label,curves=TRUE){
 V<-make_V(ff);rows<-list();ag<-list()
 add<-function(o,analysis,age=NA_real_,sex='Difference',period='2015-2018'){
  se<-sqrt(pmax(0,rowSums((o$G%*%V)*o$G)))
  rows[[length(rows)+1]]<<-data.table(model=label,analysis=analysis,age=age,sex=sex,period=period,outcome=names(o$est),estimate=o$est,se=se)
 }
 diff<-function(m,f)list(est=m$est-f$est,G=m$G-f$G)
 ref<-z[sex=='Male',.N,by=vars];ref[,weight:=N/sum(N)]
 for(pp in levels(z$period)){
  m<-predict_joint(ff,ref,'Male',pp);f<-predict_joint(ff,ref,'Female',pp)
  add(m,'overall',sex='Male',period=pp);add(f,'overall',sex='Female',period=pp);add(diff(m,f),'overall',period=pp)
 }
 if(curves){
  ref<-z[sex=='Male',.N,by=.(stage,grade,er,race)]
  ref[,`:=`(weight=N/sum(N),age=60,sex=factor('Male',levels=levels(z$sex)),period=factor('2015-2018',levels=levels(z$period)))]
  for(a in seq(40,80,by=2.5)){
   ref[,age:=a];m<-predict_joint(ff,ref,'Male');f<-predict_joint(ff,ref,'Female');g<-diff(m,f)
   add(m,'age',a,'Male');add(f,'age',a,'Female');add(g,'age',a);ag[[as.character(a)]]<-g
  }
  add(diff(ag[['75']],ag[['50']]),'age75_minus50')
  # Direct cardiovascular-minus-other-component contrasts at fixed ages.
  for(a in c(50,75)){
   g<-ag[[as.character(a)]];ks<-setdiff(causes[2:6],'cardiovascular')
   ee<-sapply(ks,function(k)g$est['cardiovascular']-g$est[k]);GG<-t(sapply(ks,function(k)g$G['cardiovascular',]-g$G[k,]))
   names(ee)<-paste0('cardiovascular_minus_',ks);add(list(est=ee,G=GG),'cause_comparison',a)
  }
 }
 rbindlist(rows)
}
