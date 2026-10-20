source('work/revision_20261019/engine.R')
args<-commandArgs(TRUE);chain<-as.integer(args[1]);maxiter<-if(length(args)>1)as.integer(args[2]) else 40
dir.create(file.path(outroot,'private_imputations'),showWarnings=FALSE)
set.seed(20261019+1000*chain);z<-load_data();targets<-c('stage','grade','er')
miss<-setNames(lapply(targets,function(v)z[[v]]=='Unknown'),targets)
for(v in targets){
 z[,(v):=droplevels(factor(as.character(get(v)),levels=setdiff(levels(z[[v]]),'Unknown')))]
 for(sx in c('Female','Male')){ii<-which(miss[[v]] & z$sex==sx);obs<-as.character(z[[v]][!miss[[v]] & z$sex==sx]);z[ii,(v):=sample(obs,length(ii),replace=TRUE)]}
}
trace<-list();checks<-list()
for(iter in seq_len(maxiter)){
 tm<-proc.time()[3];ff<-fit_models(z,draw=TRUE)
 for(v in targets){
  # Conditional covariate model uses every record, collapsed into exact duplicate
  # covariate profiles with frequency weights. No female training subsample.
  fm<-as.formula(paste(v,'~sex*(period+ns(age,knots=c(50,65,75),Boundary.knots=c(18,85))+race+',paste(setdiff(targets,v),collapse='+'),')'))
  tr<-z[,.N,by=vars]
  f<-multinom(fm,data=tr,weights=N,trace=FALSE,Hess=TRUE,MaxNWts=10000,maxit=400,decay=.001)
  stopifnot(f$convergence==0)
  b<-coef(f);if(is.null(dim(b)))b<-matrix(b,nrow=1)
  bb<-matrix(MASS::mvrnorm(1,as.vector(t(b)),vcov(f)),nrow=nrow(b),byrow=TRUE)
  ii<-which(miss[[v]]);zz<-copy(z[ii]);X<-model.matrix(delete.response(terms(f)),zz)
  eta<-cbind(0,X%*%t(bb));ll<-matrix(0,nrow(zz),length(f$lev))
  for(j in seq_along(f$lev)){
   zz[,(v):=factor(f$lev[j],levels=levels(z[[v]]))]
   ll[,j]<-loglik(ff,zz)
  }
  lp<-eta+ll;lp<-lp-apply(lp,1,max);pr<-exp(lp);pr<-pr/rowSums(pr)
  stopifnot(all(is.finite(pr)),max(abs(rowSums(pr)-1))<1e-10)
  u<-runif(nrow(pr));cp<-pr
  for(j in 2:ncol(cp))cp[,j]<-cp[,j]+cp[,j-1]
  idx<-1L+rowSums(sweep(cp[,seq_len(ncol(cp)-1),drop=FALSE],1,u,'<'))
  z[ii,(v):=f$lev[idx]]
  for(sx in c('Male','Female'))for(lev in levels(z[[v]]))trace[[length(trace)+1]]<-data.table(chain=chain,iter=iter,variable=v,sex=sx,category=lev,proportion=mean(z[[v]][miss[[v]] & z$sex==sx]==lev))
 }
 fwrite(rbindlist(trace),file.path(outroot,paste0('trace_c',chain,'.csv')))
 checks[[iter]]<-data.table(chain=chain,iter=iter,seconds=proc.time()[3]-tm,models_converged=all(sapply(ff,`[[`,'converged')))
 fwrite(rbindlist(checks),file.path(outroot,paste0('imputation_checks_c',chain,'.csv')))
 cat(format(Sys.time()),'CHAIN',chain,'ITER',iter,'SECONDS',proc.time()[3]-tm,'\n');flush.console()
 if(iter %in% seq(22,40,2)){
  label<-sprintf('smc_c%d_i%02d',chain,(iter-20)/2)
  saveRDS(z[,c('row',targets),with=FALSE],file.path(outroot,'private_imputations',paste0(label,'.rds')))
 }
}
