source('work/revision_20261019/engine.R')
d<-load_data();set.seed(2);z<-d[sample(.N,20000)]
for(v in c('stage','grade','er')){z[get(v)=='Unknown',(v):=if(v=='stage')'Localized' else if(v=='grade')'II' else 'Positive'];z[,(v):=droplevels(get(v))]}
ff<-readRDS(file.path(outroot,'fits/unknown_category.rds'))$fits
# Fit dimensions differ between unknown-category and completed covariate levels;
# use a modest random coefficient vector solely for the algebra/benchmark check.
zz<-copy(z);zz[,interval:=factor(1,levels=1:5)];p<-ncol(model.matrix(form,zz))
ff<-setNames(lapply(causes,function(k)list(b=rnorm(p,0,.04))),causes)
rows<-list()
for(v in c('stage','grade','er')){
 cats<-levels(z[[v]]);t1<-system.time({old<-sapply(cats,function(c){zz<-copy(z);zz[,(v):=factor(c,levels=cats)];loglik(ff,zz)})})[3]
 t2<-system.time(new<-category_loglik(ff,z,v,cats))[3]
 stopifnot(max(abs(old-new))<1e-10)
 rows[[v]]<-data.table(variable=v,old_seconds=t1,new_seconds=t2,max_error=max(abs(old-new)))
}
fwrite(rbindlist(rows),file.path(outroot,'category_likelihood_validation.csv'));print(rbindlist(rows))
