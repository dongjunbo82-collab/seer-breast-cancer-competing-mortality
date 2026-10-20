source('work/revision_20261019/engine.R')
z<-load_data();set.seed(1);z<-z[sample(.N,80)]
z[stage=='Unknown',stage:='Localized'];z[grade=='Unknown',grade:='II'];z[er=='Unknown',er:='Positive']
zz<-copy(z);zz[,interval:=factor(1,levels=1:5)];p<-ncol(model.matrix(form,zz))
ff<-setNames(lapply(causes,function(k)list(b=rnorm(p,0,.04))),causes)
fast<-loglik(ff,z);slow<-numeric(nrow(z))
for(i in seq_len(nrow(z)))for(j in 1:5){
 a<-copy(z[i]);a[,interval:=factor(j,levels=1:5)];eta<-sapply(ff,function(f)sum(model.matrix(form,a)*f$b));h<-exp(eta)
 if(a$period!='2015-2018' || j==1)h[5]<-0
 slow[i]<-slow[i]-sum(h)*max(0,min(a$time,j)-(j-1))
 if(a$time>j-1 && a$time<=j && a$detail %in% causes)slow[i]<-slow[i]+eta[match(a$detail,causes)]
}
stopifnot(max(abs(fast-slow))<1e-10)
# Constant hazards reduce to event log-hazard minus total cumulative exposure.
for(k in causes)ff[[k]]$b[]<-0
got<-loglik(ff,z);expected<--6*pmin(z$time,5)-ifelse(z$period=='2015-2018',pmax(0,pmin(z$time,5)-1),0)
stopifnot(max(abs(got-expected))<1e-10)
fwrite(data.table(check=c('individual_vs_vectorized','constant_hazard_closed_form'),max_error=c(max(abs(fast-slow)),max(abs(got-expected)))),file.path(outroot,'likelihood_validation.csv'))
cat('Likelihood checks passed\n')
