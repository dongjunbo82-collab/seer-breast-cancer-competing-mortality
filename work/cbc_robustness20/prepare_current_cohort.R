suppressPackageStartupMessages(library(data.table))
dir.create('work',showWarnings=FALSE)
d<-fread('work/extended.csv.gz')
d[,period:=cut(year,c(1999,2004,2009,2014,2018,2020,2023),labels=c('2000-2004','2005-2009','2010-2014','2015-2018','2019-2020','2021-2023'))]
d[,time:=pmax(months,.5)/12]
for(v in c('sex','period','race','stage','grade','er'))d[,(v):=factor(get(v))]
d[,sex:=relevel(sex,'Female')];d[,stage:=relevel(stage,'Localized')]
main<-d[mainhist==1 & year<=2018];main[,period:=droplevels(period)]
stopifnot(nrow(main)==882752,sum(main$sex=='Male')==6111)
saveRDS(main,'work/main.rds',compress=FALSE)
cat('Current cohort prepared:',nrow(main),'records. No mortality models fitted.\n')
