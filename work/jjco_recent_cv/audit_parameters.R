library(data.table);dest<-'work/jjco_final_tightening';dir.create(dest,recursive=TRUE,showWarnings=FALSE);rr<-list();tt<-list()
for(spec in c('with_sex_interval','without_sex_interval'))for(m in c('unknown_category','complete_case')){
 folder<-if(spec=='with_sex_interval')'work/jjco_recent_cv' else 'work/jjco_recent_cv/no_sex_interval';x<-readRDS(file.path(folder,paste0(m,'_private_fit.rds')))
 for(i in seq_along(x$fits)){
  f<-x$fits[[i]];b<-f$b;se<-sqrt(diag(f$V));cv<-cov2cor(f$V);diag(cv)<-0;cause<-c('breast','cardiovascular','remaining_deaths')[i];male<-grepl('sexMale',names(b));it<-grepl('sexMale:interval|interval.*:sexMale',names(b));nonint<-names(b)!='(Intercept)'
  rr[[length(rr)+1]]<-data.table(specification=spec,model=m,cause=cause,n_parameters=length(b),male_terms=sum(male),male_age_terms=sum(grepl('sexMale:age10',names(b))),male_interval_terms=sum(it),max_abs_nonintercept=max(abs(b[nonint])),max_se=max(se),max_abs_correlation=max(abs(cv)),min_covariance_eigenvalue=min(eigen(f$V,symmetric=TRUE,only.values=TRUE)$values),finite_parameters=all(is.finite(b)),finite_se=all(is.finite(se)))
  tt[[length(tt)+1]]<-data.table(specification=spec,model=m,cause=cause,term=names(b),estimate=as.numeric(b),se=se,lower=b-1.96*se,upper=b+1.96*se,male_related=male)
 }
}
fwrite(rbindlist(rr),file.path(dest,'parameter_audit.csv'));fwrite(rbindlist(tt),file.path(dest,'coefficient_audit.csv'));print(rbindlist(rr)[cause=='cardiovascular'])
