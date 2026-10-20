import os
import pathlib,csv,re,zipfile,xml.etree.ElementTree as E,xlrd
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def xlsxrows(f):
 z=zipfile.ZipFile(f);ss=[''.join(n.itertext()) for n in E.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
 for r in E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:row',ns):
  cells={}
  for c in r.findall('s:c',ns):
   v=c.find('s:v',ns);t=v.text if v is not None else '';cells[re.sub(r'\d','',c.get('r'))]=ss[int(t)] if c.get('t')=='s' and t else t
  yield [cells.get(k,'') for k in ['A','B','C','D']]
result=[]
for r in csv.DictReader(open('work/lifetables/manifest.csv',encoding='utf8')):
 f=r['file'];y=int(r['year']);sex=r['sex']
 if f.endswith('.xls'):
  s=xlrd.open_workbook(f).sheet_by_index(0);rows=[s.row_values(i) for i in range(s.nrows)]
 else:rows=list(xlsxrows(f))
 title=' '.join(str(x) for row in rows[:4] for x in row)
 assert str(y) in title and (('females' in title.lower()) if sex=='Female' else ('males' in title.lower() and 'females' not in title.lower())),(f,title)
 found={}
 for row in rows:
  age=re.match(r'^\s*(\d{1,3})(?:\D|$)',str(row[0]))
  if not age:continue
  try:q=float(row[1])
  except (ValueError,TypeError):continue
  if not 0<=q<=1:continue
  a=int(age[1]);assert a not in found,(f,a);found[a]=q
 assert all(a in found for a in range(100)),(f,len(found))
 result.extend(dict(year=y,sex=sex,age=a,q=q) for a,q in found.items())
with open('work/lifetables/qx.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['year','sex','age','q']);w.writeheader();w.writerows(result)
print('validated',len(result),'life-table cells')
