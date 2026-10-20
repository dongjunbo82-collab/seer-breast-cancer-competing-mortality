source('work/revision_20261019/engine.R')
files<-list.files(file.path(outroot,'private_imputations'),pattern='^smc_.*[.]rds$',full.names=TRUE)
stopifnot(length(files)==40);d<-load_data();rows<-list()
for(fn in files){
 a<-readRDS(fn);stopifnot(identical(a$row,d$row))
 for(v in c('stage','grade','er')){known<-d[[v]]!='Unknown';stopifnot(all(as.character(a[[v]][known])==as.character(d[[v]][known])),!anyNA(a[[v]]),!any(as.character(a[[v]])=='Unknown'))}
 label<-sub('[.]rds$','',basename(fn));o<-readRDS(file.path(outroot,'fits',paste0(label,'.rds')))
 for(k in causes){f<-o$fits[[k]];rows[[length(rows)+1]]<-data.table(imputation=label,cause=k,events=f$events,max_abs_coefficient=max(abs(f$b)),max_SE=max(sqrt(diag(f$V))))}
}
fwrite(rbindlist(rows),file.path(outroot,'coefficient_diagnostics.csv'))
print(rbindlist(rows)[,.(max_coef=max(max_abs_coefficient),max_SE=max(max_SE)),by=cause])
cat('40 imputations: row alignment, observed entries and categories verified\n')
