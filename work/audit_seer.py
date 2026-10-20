import os
import csv,json,collections,pathlib
src=pathlib.Path(os.environ["SEER_BASE_EXPORT"])
out=pathlib.Path(__file__).resolve().parent;out.mkdir(exist_ok=True)
counts=collections.Counter(); groups=collections.defaultdict(collections.Counter); values=collections.defaultdict(collections.Counter)
keys=set(); dup=0; bad=0
with src.open(encoding='utf-8-sig',newline='') as f:
 r=csv.reader(f); header=next(r)
 for row in r:
  if len(row)!=len(header): bad+=1;continue
  sex=row[3]; year=int(row[4]); period=f'{2000+(year-2000)//5*5}-{min(2023,2004+(year-2000)//5*5)}'
  counts[sex]+=1;g=groups[(sex,period)];g['records']+=1
  key=(row[0],row[2]);dup+=key in keys;keys.add(key)
  for i in range(3,30):
   if i not in [4,5,6,10,11,23]:values[header[i]][row[i]]+=1
  dead=row[25]=='Dead';bc=row[26]=='Dead (attributable to this cancer dx)';other=row[27]=='Dead (attributable to causes other than this cancer dx)'
  g['dead']+=dead;g['bc_death']+=bc;g['other_death']+=other;g['dead_unclassified']+=dead and not(bc or other)
  try:t=int(row[23])
  except ValueError:t=9999
  g['survival_unknown']+=t==9999;g['zero_months']+=t==0
  g['bc_death_before60']+=bc and t<60;g['other_death_before60']+=other and t<60
  g['alive_censored_before60']+=(not dead) and t<60
  g['sequence_first_or_only']+=row[1] in ['One primary only','1st of 2 or more primaries']
  stage=row[15] if year<2004 else row[16]
  g['stage_unknown']+=('unknown' in stage.lower() or stage in ['Blank(s)',''])
  grade=row[17] if year<2018 else row[18]
  g['grade_unknown']+=('unknown' in grade.lower() or grade in ['Blank(s)',''])
  g['er_unknown']+=row[19] not in ['Positive','Negative']
  if year>=2010:g['subtype_unknown']+=row[22] not in ['HR+/HER2-','HR+/HER2+','HR-/HER2+','HR-/HER2-']
result=dict(source=str(src),columns=header,counts=dict(counts),malformed_rows=bad,duplicate_patient_record_keys=dup,groups=[dict(sex=k[0],period=k[1],**v) for k,v in sorted(groups.items())],values=values)
(out/'seer_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ['values','columns']},indent=2))
print(json.dumps({k:v for k,v in values.items() if k in [header[25],header[26],header[27],header[28],header[18],header[22]]},indent=2))
