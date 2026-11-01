# Set SEER_ANALYSIS_ROOT to the locally authorized input project.
source('work/revision_20261019/engine.R')
ee<-new.env(parent=globalenv());ex<-parse(file.path(source_root,'work/models.R'));for(e in ex)if(is.call(e)&&as.character(e[[1]])%in%c('<-','=')&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name('function')))eval(e,ee)
setDTthreads(1);dest<-'work/jjco_robustness_final';d<-load_data();targets<-c('stage','grade','er');rr<-list();checks<-list()
# Three causes use a shared age spline, deliberately simplifying the original model.
fm<-~sex+ns(age,knots=c(50,65,75),Boundary.knots=c(18,85))+interval+sex:interval+stage+stage:interval+grade+er+race
calc<-function(z,label,exclude60=FALSE){
 if(exclude60)z[months==60,event:=0L]
 z<-droplevels(z);ff<-ee$fit_hazards(z,fm);ref<-z[sex=='Male',.N,by=vars];ref[,weight:=N/sum(N)];V<-matrix(0,3*length(ff[[1]]$b),3*length(ff[[1]]$b));p<-length(ff[[1]]$b)
 for(k in 1:3){ix<-((k-1)*p+1):(k*p);V[ix,ix]<-ff[[k]]$V}
 m<-ee$predict_risk(ff,ref,'Male',levels(z$period)[1]);f<-ee$predict_risk(ff,ref,'Female',levels(z$period)[1])
 ans<-list(Male=m,Female=f,Difference=list(est=m$est-f$est,G=m$G-f$G))
 for(ss in names(ans)){a<-ans[[ss]];se<-sqrt(pmax(0,rowSums((a$G%*%V)*a$G)));rr[[length(rr)+1]]<<-data.table(model=label,sex=ss,outcome=names(a$est),estimate=a$est,se=se,lower=a$est-1.96*se,upper=a$est+1.96*se)}
 checks[[length(checks)+1]]<<-data.table(model=label,n_male=sum(z$sex=='Male'),n_female=sum(z$sex=='Female'),male_deaths_before60=sum(z$sex=='Male'&z$event>0&z$months<60),male_deaths_at60=sum(z$sex=='Male'&z$event>0&z$months==60),converged=all(sapply(ff,`[[`,'converged')))
 fwrite(rbindlist(rr),file.path(dest,'results.csv'));fwrite(rbindlist(checks),file.path(dest,'counts.csv'));cat('DONE',label,'\n');flush.console()
}
for(mode in c('unknown_category','complete_case')){
 z<-copy(d[year>=2015]);if(mode=='complete_case')z<-z[stage!='Unknown'&grade!='Unknown'&er!='Unknown'];calc(copy(z),paste0('recent_',mode));calc(copy(z),paste0('recent_exclude60_',mode),TRUE)
 z<-copy(d[year>=2010&year<=2014]);if(mode=='complete_case')z<-z[stage!='Unknown'&grade!='Unknown'&er!='Unknown'];calc(z,paste0('prepandemic5_',mode),TRUE)
}
fwrite(d[,.(n=.N,stage_missing=sum(stage=='Unknown'),grade_missing=sum(grade=='Unknown'),er_missing=sum(er=='Unknown')),by=.(sex,period)],file.path(dest,'missingness.csv'))
