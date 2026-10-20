Sys.setenv(LOAD_ONLY='1');source('work/models.R')
setDTthreads(3)
# Common three-year horizon ending before 2020 for every eligible diagnosis year.
z<-copy(main[year<=2016]);z[,period:=factor(cut(year,c(1999,2004,2009,2014,2016),labels=c('2000-2004','2005-2009','2010-2014','2015-2016')))]
# Events recorded at month 36 are censored at the boundary: monthly rounding
# must not carry an event into 2020 for a late-December 2016 diagnosis.
z[months>=36,`:=`(time=3,event=0L)]
analyse(z,'revision_prepandemic3',horizon=3)
analyse(z[stage!='Unknown' & grade!='Unknown' & er!='Unknown'],'revision_prepandemic3_complete',horizon=3)
# Independently fitted age-stratum models permit all coefficients to differ.
young<-copy(main[age<65]);old<-copy(main[age>=65])
fy<-as.formula("~sex*period+ns(age,knots=c(40,55),Boundary.knots=c(18,65))+interval+sex:interval+stage+stage:interval+grade+er+race")
fo<-as.formula("~sex*period+ns(age,knots=c(70,78),Boundary.knots=c(65,85))+interval+sex:interval+stage+stage:interval+grade+er+race")
analyse(young,'revision_age_under65',fy)
analyse(old,'revision_age_65plus',fo)
labels<-c('revision_prepandemic3','revision_prepandemic3_complete','revision_age_under65','revision_age_65plus')
rr<-rbindlist(lapply(labels,function(x)fread(paste0('work/models/',x,'_results.csv'))))
fwrite(rr,'outputs/tables/revision_additional_results.csv')
# Independent patients across age strata; fixed, stratum-specific references.
a<-rr[model=='revision_age_under65' & kind=='gap'];b<-rr[model=='revision_age_65plus' & kind=='gap']
cmp<-merge(a,b,by=c('period','outcome'),suffixes=c('_young','_older'))
cmp[,`:=`(difference=estimate_older-estimate_young,contrast_se=sqrt(se_older^2+se_young^2))]
cmp[,`:=`(lower=difference-1.96*contrast_se,upper=difference+1.96*contrast_se,p=2*pnorm(-abs(difference/contrast_se)))]
cmp[,p_holm:=p.adjust(p,'holm')]
fwrite(cmp[,.(period,outcome,difference,contrast_se,lower,upper,p,p_holm)],'outputs/tables/revision_age_contrasts.csv')
counts<-rbindlist(lapply(list(under65=young,age65plus=old,prepandemic3=z),function(x){h<-if(identical(x,z))3 else 5;x[,.(N=.N,bc=sum(event==1 & time<=h),other=sum(event==2 & time<=h),unknown=sum(event==3 & time<=h)),by=.(sex,period)]}),idcol='analysis')
fwrite(counts,'outputs/tables/revision_additional_counts.csv')
stopifnot(all(is.finite(rr$estimate)),all(rr$lower<=rr$upper))
cat('REVISION ANALYSES COMPLETE\n')
