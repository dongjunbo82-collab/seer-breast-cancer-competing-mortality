source('work/revision_20261019/engine.R')
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
dir.create(file.path(outroot,'figures'),showWarnings=FALSE)
labs<-c(all='All-cause',breast='Breast cancer',other_total='Other causes combined',cardiovascular='Cardiovascular',other_cancer='Other neoplasms',respiratory='COPD or pneumonia/influenza',covid='COVID-19',remaining_other='Remaining other causes',unknown='Unknown cause')
savefig<-function(g,name,w,h){
 fn<-file.path(outroot,'figures',name)
 svglite::svglite(paste0(fn,'.svg'),width=w,height=h);print(g);dev.off()
 ragg::agg_tiff(paste0(fn,'.tif'),width=w,height=h,units='in',res=600,compression='lzw');print(g);dev.off()
 ragg::agg_png(paste0(fn,'.png'),width=w,height=h,units='in',res=180);print(g);dev.off()
}
theme_set(theme_classic(base_size=10,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold',size=10),legend.position='bottom',legend.title=element_blank(),plot.caption=element_text(hjust=0,size=8),plot.margin=margin(7,10,7,7)))
x<-fread(file.path(outroot,'revised_MI.csv'));ss<-fread(file.path(outroot,'revised_sensitivity.csv'))
oo<-c('all','breast','other_total','cardiovascular','other_cancer','respiratory','covid','remaining_other','unknown')
q<-x[analysis=='overall' & sex=='Difference' & period=='2015-2018'];q[,label:=factor(labs[outcome],levels=rev(labs[oo]))]
q[,type:=ifelse(outcome %in% c('all','other_total'),'Aggregate','Component')]
g<-ggplot(q,aes(100*estimate,label,color=type))+geom_vline(xintercept=0,color='grey65',linewidth=.4)+geom_errorbar(aes(xmin=100*lower,xmax=100*upper),orientation='y',width=.18,linewidth=.6)+geom_point(size=2.3)+scale_color_manual(values=c(Aggregate='black',Component='#287FA3'))+labs(x='Male minus female five-year risk (percentage points)',y=NULL,caption='2015-2018; common pooled male case mix; 40 imputations. Bars: 95% confidence intervals.\nOther causes combined is the sum of its five detailed components; aggregate rows are not additional causes.')+theme(legend.position='none')
savefig(g,'Fig1',7.4,4.5)
q<-x[analysis=='age' & outcome %in% c('all','breast','other_total')];q[,label:=factor(labs[outcome],levels=labs[c('all','breast','other_total')])]
q1<-q[sex!='Difference'];q2<-q[sex=='Difference']
a<-ggplot(q1,aes(age,100*estimate,color=sex,fill=sex))+geom_ribbon(aes(ymin=100*lower,ymax=100*upper),alpha=.14,color=NA)+geom_line(linewidth=.8)+facet_wrap(~label,nrow=1)+scale_color_manual(values=c(Female='#C87730',Male='#287FA3'))+scale_fill_manual(values=c(Female='#C87730',Male='#287FA3'))+labs(x=NULL,y='Five-year mortality risk (%)')
b<-ggplot(q2,aes(age,100*estimate))+geom_hline(yintercept=0,color='grey65',linewidth=.4)+geom_ribbon(aes(ymin=100*lower,ymax=100*upper),fill='#287FA3',alpha=.16)+geom_line(color='#287FA3',linewidth=.8)+facet_wrap(~label,nrow=1)+labs(x='Age at diagnosis (years)',y='Risk difference (percentage points)')
g<-(a/b)+plot_annotation(tag_levels='A',caption='2015-2018 predictions at one fixed measured non-age male case mix. Shading: pointwise 95% intervals.\nAge curves are exploratory; sparse covariate combinations at the age extremes require caution.')
savefig(g,'Fig2',7.4,6.1)
q<-x[analysis=='age' & sex=='Difference' & outcome %in% c('cardiovascular','other_cancer','respiratory','remaining_other')]
q[,label:=factor(labs[outcome],levels=labs[c('cardiovascular','other_cancer','respiratory','remaining_other')])]
g<-ggplot(q,aes(age,100*estimate))+geom_hline(yintercept=0,color='grey65',linewidth=.4)+geom_ribbon(aes(ymin=100*lower,ymax=100*upper),fill='#287FA3',alpha=.16)+geom_line(color='#287FA3',linewidth=.8)+facet_wrap(~label,ncol=2,scales='free_y')+labs(x='Age at diagnosis (years)',y='Risk difference (percentage points)',caption='Male minus female five-year risk; 2015-2018 at a fixed measured non-age case mix.\nShading: pointwise 95% intervals. Panel scales differ. COVID-19 and unknown-cause estimates are in the tables;\nonly 11 and 24 male deaths, respectively, precluded estimating sex-specific age splines for these causes.')
savefig(g,'Fig3',7.4,5.8)
q<-rbind(x[analysis=='age' & sex=='Difference' & outcome %in% c('all','breast','other_total'),.(age,outcome,estimate,model='Likelihood-weighted FCS')],ss[analysis=='age' & sex=='Difference' & outcome %in% c('all','breast','other_total'),.(age,outcome,estimate,model=ifelse(model=='complete_case','Complete cases','Unknown categories'))])
q[,label:=factor(labs[outcome],levels=labs[c('all','breast','other_total')])]
g<-ggplot(q,aes(age,100*estimate,color=model,linetype=model))+geom_hline(yintercept=0,color='grey65')+geom_line(linewidth=.75)+facet_wrap(~label,nrow=1)+scale_color_manual(values=c('Likelihood-weighted FCS'='#287FA3','Complete cases'='#C87730','Unknown categories'='#8064A2'))+labs(x='Age at diagnosis (years)',y='Risk difference (percentage points)',caption='Each analysis fixes its own non-age reference across ages. Complete cases change the eligible reference population.\nAll three specifications use the same seven-cause outcome models; intervals and fixed contrasts are tabulated separately.')
savefig(g,'S13',7.4,3.6)
tr<-rbindlist(lapply(1:4,function(k)fread(file.path(outroot,paste0('trace_c',k,'.csv')))))
tr<-tr[(variable=='stage' & category=='Regional') | (variable=='grade' & category=='III') | (variable=='er' & category=='Positive')]
tr[,label:=paste(variable,category)];tr[,chain:=factor(chain)]
g<-ggplot(tr,aes(iter,proportion,color=chain))+geom_vline(xintercept=20.5,color='grey60',linetype=2)+geom_line(linewidth=.5)+facet_grid(label~sex,scales='free_y')+scale_color_manual(values=c('#287FA3','#C87730','#8064A2','#4D8556'))+labs(x='Sweep',y='Proportion among imputed entries',caption='Four independent chains. Dashed line marks the end of burn-in; sweeps 22, 24, ..., 40 were retained.\nAll category diagnostics, rather than only these representative traces, are reported in Table S29.')
savefig(g,'S16',7.4,6)
su<-fread(file.path(outroot,'support_grid.csv'))
q<-melt(su,id.vars=c('scope','age'),measure.vars=c('zero_weight','below5_weight'),variable.name='criterion',value.name='weight')
q[,criterion:=factor(criterion,levels=c('below5_weight','zero_weight'),labels=c('Fewer than 5 cases in either sex','Zero cases in either sex'))]
g<-ggplot(q,aes(age,100*weight,color=criterion,linetype=criterion))+geom_line(linewidth=.8)+facet_wrap(~scope,nrow=1)+scale_color_manual(values=c('#C87730','#287FA3'))+labs(x='Prediction age (years)',y='Reference weight (%)',caption='Support assessed in the 10-year window centered on each prediction age, using stage, grade, ER, and race profiles.\nValues average across 40 completed datasets; support is descriptive and does not establish model validity.')
savefig(g,'S17',7.4,3.6)
