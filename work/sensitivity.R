Sys.setenv(LOAD_ONLY='1');source('work/models.R')
analyse_new<-analyse
analyse<-function(z,label,...) {if(file.exists(paste0('work/models/',label,'_results.csv')))return(invisible(NULL));analyse_new(z,label,...)}
analyse(main[age<85],'age_under85')
analyse(main[year<=2017],'before_grade_change')
analyse(main[stage!='Unknown' & grade!='Unknown' & er!='Unknown'],'complete_case')
analyse(d[year<=2018],'all_histology')
analyse(main[micro==1],'microscopic')
analyse(main[months>0],'exclude_zero')
analyse(main[stage %in% c('Localized','Regional')],'nonmetastatic')
analyse(main,'sex_specific',as.formula(paste(deparse(forms$full,width.cutoff=500),'+ sex:stage + sex:ns(age,knots=c(50,65,75),Boundary.knots=c(18,85))')))
vars0<-vars;vars<-c(vars,'pr');main[,pr:=factor(pr)]
analyse(main,'pr_adjusted',as.formula(paste(deparse(forms$full,width.cutoff=500),'+ pr')))
vars<-vars0
z<-d[mainhist==1 & year<=2020];analyse(z,'three_year',horizon=3)
z<-main[year<=2014];analyse(z,'prepandemic_followup')
z<-main[year>=2004];analyse(z,'stage_2004on')
z<-main[year>=2010 & subtype=='HR+/HER2-'];analyse(z,'hrpos_her2neg')
cat('ALL SENSITIVITIES COMPLETE\n')
