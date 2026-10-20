Sys.setenv(LOAD_ONLY='1');source('work/models.R')
q<-fread('work/lifetables/qx.csv');setkey(q,sex,year,age)
ref<-main[sex=='Male' & age<85,.N,by=age];ref[,w:=N/sum(N)]
results<-list()
for(pp in levels(main$period))for(ss in c('Male','Female')){
 years<-sort(unique(main[period==pp]$year));a<-CJ(age=ref$age,year=years);a<-merge(a,ref[,.(age,w)],by='age');a[,w:=w/length(years)]
 for(offset in c(0,.5)){
  surv<-rep(1,nrow(a))
  # Cohort progression: age and calendar year increase annually. The half-year
  # alternative averages adjacent calendar-year hazards; no fabricated birthday.
  for(j in 0:4){qx<-q[data.table(sex=ss,year=a$year+j,age=a$age+j)]$q
   qnext<-q[data.table(sex=ss,year=a$year+j+1,age=a$age+j)]$q
   surv<-surv*(1-qx)^(1-offset)*(1-qnext)^offset
  }
  stopifnot(!anyNA(surv));results[[length(results)+1]]<-data.table(period=pp,sex=ss,offset=offset,expected_death=sum(a$w*(1-surv)))
 }
}
bg<-rbindlist(results);fwrite(bg,'outputs/tables/population_expected_mortality.csv')
gap<-dcast(bg,period+offset~sex,value.var='expected_death');gap[,expected_gap:=Male-Female]
obs<-fread('work/models/age_under85_results.csv')[kind=='gap' & outcome=='all',.(period,observed_gap=estimate,observed_lower=lower,observed_upper=upper)]
gap<-merge(gap,obs,by='period');gap[,descriptive_difference:=observed_gap-expected_gap]
fwrite(gap,'outputs/tables/population_benchmark.csv')
cat('Population benchmark completed\n');print(gap)
