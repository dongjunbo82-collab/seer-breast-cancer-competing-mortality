from pathlib import Path
import csv,gzip,json,statistics,collections
root=Path(__file__).resolve().parent.parent
values=collections.defaultdict(list)
with gzip.open(root/'work/extended.csv.gz','rt',encoding='utf-8') as f:
 for row in csv.DictReader(f):
  if row['mainhist'] != '1' or int(row['year']) > 2018: continue
  values[row['sex']].append(min(60,max(0.5,float(row['months'])))/12)
summary={sex:{'n':len(v),'total_person_years':sum(v),'mean_years':statistics.mean(v),'median_years':statistics.median(v)} for sex,v in values.items()}
assert summary['Male']['n']==6111 and summary['Female']['n']==876641
(root/'outputs/tables/observed_followup_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
