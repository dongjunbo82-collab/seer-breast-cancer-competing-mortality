import os
import csv,gzip,pathlib,collections,json,hashlib
root=pathlib.Path(__file__).resolve().parent.parent
src=pathlib.Path(os.environ["SEER_BASE_EXPORT"])
counts=collections.Counter()
with src.open(encoding='utf-8-sig',newline='') as f,gzip.open(root/'work/extended.csv.gz','wt',encoding='utf8',newline='') as o:
 rd=csv.DictReader(f);wr=csv.writer(o)
 wr.writerow(['row','sex','year','age','ageband','race','stage','grade','er','pr','her2','subtype','months','event','hist','first','micro','mainhist','rawgrade','marital'])
 for line,r in enumerate(rd,2):
  if r['Type of Reporting Source'] in ['Death certificate only','Autopsy only']:continue
  try:t=int(r['Survival months'])
  except:continue
  if t>=9999:continue
  y=int(r['Year of diagnosis']);h=int(r['Histologic Type ICD-O-3'])
  first=r['Sequence number'] in ['One primary only','1st of 2 or more primaries']
  if not first:continue
  st=r['Summary stage 2000 (1998-2017)'] if y<2004 else r['Combined Summary Stage with Expanded Regional Codes (2004+)']
  stage=next((x for x in ['Localized','Regional','Distant'] if st.startswith(x)),'Unknown')
  raw=r['Grade Recode (thru 2017)'] if y<2018 else r['Derived Summary Grade 2018 (2018+)']
  gm={'Well differentiated; Grade I':'I','Moderately differentiated; Grade II':'II','Poorly differentiated; Grade III':'III','Undifferentiated; anaplastic; Grade IV':'III','Site-specific grade system category (1)':'I','Site-specific grade system category (2)':'II','Site-specific grade system category (3)':'III','Well differentiated':'I','Moderately differentiated':'II','Poorly differentiated':'III','Undifferentiated and anaplastic':'III'}
  grade=gm.get(raw,'Unknown')
  ev=0 if r['Vital status recode (study cutoff used)']=='Alive' else 1 if r['SEER cause-specific death classification']=='Dead (attributable to this cancer dx)' else 2 if r['SEER other cause of death classification']=='Dead (attributable to causes other than this cancer dx)' else 3
  er=r['ER Status Recode Breast Cancer (1990+)'];pr=r['PR Status Recode Breast Cancer (1990+)']
  if er not in ['Positive','Negative']:er='Unknown'
  if pr not in ['Positive','Negative']:pr='Unknown'
  race=r['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']
  race={'Non-Hispanic White':'NHW','Non-Hispanic Black':'NHB','Hispanic (All Races)':'Hispanic','Non-Hispanic Asian or Pacific Islander':'NHAPI','Non-Hispanic American Indian/Alaska Native':'Other'}.get(race,'Other')
  marital=r['Marital status at diagnosis']
  marital='Married' if marital=='Married (including common law)' else 'Widowed' if marital=='Widowed' else 'Single' if marital=='Single (never married)' else 'Separated_divorced' if marital in ['Separated','Divorced'] else 'Unknown'
  wr.writerow([line,r['Sex'],y,int(r['Age recode with single ages and 85+'].split()[0].replace('+','')),r['Age recode with <1 year olds and 90+'],race,stage,grade,er,pr,r['Derived HER2 Recode (2010+)'],r['Breast Subtype (2010+)'],t,ev,h,1,int(r['Diagnostic Confirmation']=='Positive histology'),int(8010<=h<=8589),raw,marital])
  counts[(r['Sex'],int(8010<=h<=8589),y<=2018)]+=1
print(counts)
