source('work/revision_20261019/engine.R')
dir.create(Sys.getenv('SEER_EXTENSION_OUT','outputs/chain_stability'),recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(data.table))
d<-load_data();imp<-readRDS(file.path(outroot,'private_imputations','smc_c1_i10.rds'));for(v in c('stage','grade','er'))d[,(v):=imp[[v]]]
ff<-readRDS(file.path(outroot,'fits','smc_c1_i10.rds'))$fits;set.seed(20261020);small<-d[sample(.N,200)];rows<-list()
for(v in c('stage','grade','er')){
 cats<-levels(d[[v]]);a<-category_loglik(ff,small,v,cats);b<-sapply(cats,function(lev){z<-copy(small);z[,(v):=factor(lev,levels=levels(d[[v]]))];loglik(ff,z)});err<-max(abs(a-b));stopifnot(err<1e-10);rows[[v]]<-data.table(variable=v,n=200,maximum_absolute_loglik_difference=err)
}
fwrite(rbindlist(rows),file.path(Sys.getenv('SEER_EXTENSION_OUT','outputs/chain_stability'),'continuation_kernel_validation.csv'));print(rbindlist(rows))
