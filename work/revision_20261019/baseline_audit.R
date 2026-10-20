source('work/revision_20261019/engine.R')
rr<-list()
for(label in c('unknown_category','complete_case')){
 o<-readRDS(file.path(outroot,'fits',paste0(label,'.rds')))
 for(k in causes){f<-o$fits[[k]];rr[[length(rr)+1]]<-data.table(model=label,cause=k,max_coef=max(abs(f$b)),max_SE=max(sqrt(diag(f$V))))}
}
fwrite(rbindlist(rr),file.path(outroot,'baseline_coefficient_diagnostics.csv'));print(rbindlist(rr))
