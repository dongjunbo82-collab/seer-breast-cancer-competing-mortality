import os
import urllib.request, pathlib, re, concurrent.futures, csv, time
def fetch(url):
 for attempt in range(4):
  try:return urllib.request.urlopen(url,timeout=45).read()
  except Exception:
   if attempt==3:raise
   time.sleep(1)
out=pathlib.Path('work/lifetables');out.mkdir(exist_ok=True)
vol={2000:'58_10',2001:'52_14',2002:'53_06',2003:'54_14',2004:'56_09',2005:'58_10',2006:'58_21',2007:'59_09',2008:'61_03',2009:'62_07',2010:'63_07',2011:'64_11',2012:'65_08',2013:'66_03',2014:'66_04',2015:'67_07',2016:'68_04',2017:'68_07',2018:'69-12',2019:'70-19',2020:'71-01',2021:'72-12',2022:'74-02',2023:'74-06'}
base='https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/'
def run(item):
 y,v=item; s=fetch(base+v+'/').decode();links=re.findall('HREF="([^"]+)"',s);result=[]
 for sex,num in [('Male',2),('Female',3)]:
  if y==2000: candidates=[p for p in links if p.endswith('2000Table'+('II' if num==2 else 'III')+'.xls')]
  else:
   candidates=[p for p in links if re.search(r'table0?'+str(num)+r'(?:_Intercensal)?\.xlsx?$',p,re.I)]
   candidates.sort(key=lambda p: ('intercensal' in p.lower(),str(y) in p),reverse=True)
  if not candidates:raise ValueError((y,links))
  url='https://ftp.cdc.gov'+candidates[0];dest=out/(f'{y}_{sex}'+pathlib.Path(candidates[0]).suffix)
  if not dest.exists():dest.write_bytes(fetch(url))
  result.append({'year':y,'sex':sex,'file':str(dest),'url':url})
 print('downloaded',y,flush=True);return result
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
 rows=[r for result in ex.map(run,vol.items()) for r in result]
with (out/'manifest.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=['year','sex','file','url']);w.writeheader();w.writerows(rows)
