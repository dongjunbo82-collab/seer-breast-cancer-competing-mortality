source('work/revision_20261019/engine.R')
args<-commandArgs(TRUE);chain<-as.integer(args[1]);end<-if(length(args)>1)as.integer(args[2])else 120L
setDTthreads(1);dest<-Sys.getenv('SEER_EXTENSION_OUT','outputs/chain_stability');dir.create(file.path(dest,'private_checkpoints'),recursive=TRUE,showWarnings=FALSE)
z<-load_data();targets<-c('stage','grade','er');miss<-setNames(lapply(targets,function(v)z[[v]]=='Unknown'),targets)
statefile<-file.path(dest,'private_checkpoints',paste0('c',chain,'.rds'));tracefile<-file.path(dest,paste0('trace_c',chain,'.csv'));resultfile<-file.path(dest,paste0('estimates_c',chain,'.csv'));timingfile<-file.path(dest,paste0('timing_c',chain,'.csv'))
if(file.exists(statefile)){
 st<-readRDS(statefile);for(v in targets)z[,(v):=st$values[[v]]];assign('.Random.seed',st$rng,envir=.GlobalEnv);start<-st$iter
}else{
 st<-readRDS(file.path(outroot,'private_imputations',sprintf('smc_c%d_i10.rds',chain)));stopifnot(identical(st$row,z$row));for(v in targets)z[,(v):=st[[v]]];set.seed(20261020L+chain*100000L);start<-40L
}
traces<-if(file.exists(tracefile))list(fread(tracefile)[iter<=start])else list();res<-if(file.exists(resultfile))list(fread(resultfile)[iter<=start])else list();times<-if(file.exists(timingfile))list(fread(timingfile)[iter<=start])else list()
key_results<-function(ff,z){
 V<-make_V(ff);rr<-list();add<-function(o,analysis,age=NA_real_){rr[[length(rr)+1]]<<-data.table(analysis=analysis,age=age,outcome=names(o$est),estimate=o$est,se=sqrt(pmax(0,rowSums((o$G%*%V)*o$G))))};minus<-function(m,f)list(est=m$est-f$est,G=m$G-f$G)
 ref<-z[sex=='Male',.N,by=vars];ref[,weight:=N/sum(N)];add(minus(predict_joint(ff,ref,'Male'),predict_joint(ff,ref,'Female')),'overall')
 ref<-z[sex=='Male',.N,by=.(stage,grade,er,race)];ref[,`:=`(weight=N/sum(N),age=50,sex=factor('Male',levels=levels(z$sex)),period=factor('2015-2018',levels=levels(z$period)))];gg<-list()
 for(a in c(50,75)){ref[,age:=a];g<-minus(predict_joint(ff,ref,'Male'),predict_joint(ff,ref,'Female'));gg[[as.character(a)]]<-g;add(g,'age',a)}
 add(minus(gg[['75']],gg[['50']]),'age75_minus50');rbindlist(rr)
}
if(end<=start)quit(status=0)
for(iter in seq.int(start+1L,end)){
 tm<-proc.time()[3];ff<-fit_models(z,draw=TRUE)
 for(v in targets){
  fm<-as.formula(paste(v,'~sex*(period+ns(age,knots=c(50,65,75),Boundary.knots=c(18,85))+race+',paste(setdiff(targets,v),collapse='+'),')'));tr<-z[,.N,by=vars]
  f<-multinom(fm,data=tr,weights=N,trace=FALSE,Hess=TRUE,MaxNWts=10000,maxit=400,decay=.001);stopifnot(f$convergence==0)
  b<-coef(f);if(is.null(dim(b)))b<-matrix(b,nrow=1);bb<-matrix(MASS::mvrnorm(1,as.vector(t(b)),vcov(f)),nrow=nrow(b),byrow=TRUE)
  ii<-which(miss[[v]]);zz<-copy(z[ii]);X<-model.matrix(delete.response(terms(f)),zz);eta<-cbind(0,X%*%t(bb));ll<-category_loglik(ff,zz,v,f$lev)
  lp<-eta+ll;lp<-lp-apply(lp,1,max);pr<-exp(lp);pr<-pr/rowSums(pr);stopifnot(all(is.finite(pr)),max(abs(rowSums(pr)-1))<1e-10)
  u<-runif(nrow(pr));cp<-pr;for(j in 2:ncol(cp))cp[,j]<-cp[,j]+cp[,j-1];idx<-1L+rowSums(sweep(cp[,seq_len(ncol(cp)-1),drop=FALSE],1,u,'<'));z[ii,(v):=f$lev[idx]]
  for(sx in c('Male','Female'))for(lev in levels(z[[v]]))traces[[length(traces)+1]]<-data.table(chain=chain,iter=iter,variable=v,sex=sx,category=lev,proportion=mean(z[[v]][miss[[v]] & z$sex==sx]==lev))
 }
 if(iter>=62 && iter%%2==0){fe<-fit_models(z,draw=FALSE);ans<-key_results(fe,z);ans[,`:=`(chain=chain,iter=iter)];res[[length(res)+1]]<-ans;fwrite(rbindlist(res),resultfile)}
 times[[length(times)+1]]<-data.table(chain=chain,iter=iter,seconds=proc.time()[3]-tm)
 fwrite(rbindlist(traces),tracefile);fwrite(rbindlist(times),timingfile)
 saveRDS(list(values=z[,c('row',targets),with=FALSE],rng=.Random.seed,iter=iter),statefile)
 cat(format(Sys.time()),'CHAIN',chain,'ITER',iter,'SECONDS',proc.time()[3]-tm,'\n');flush.console()
}
