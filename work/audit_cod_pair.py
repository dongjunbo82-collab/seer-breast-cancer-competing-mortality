import os
import csv,json,itertools
from pathlib import Path
old=Path(os.environ["SEER_BASE_EXPORT"])
new=Path(os.environ["SEER_COD_EXPORT"])
with old.open(encoding='utf-8-sig',newline='') as a,new.open(encoding='utf-8-sig',newline='') as b:
    aa=csv.reader(a);bb=csv.reader(b);ha=next(aa);hb=next(bb);assert ha==hb[:30]
    n=0;changes=0
    for x,y in itertools.zip_longest(aa,bb):
        assert x is not None and y is not None
        assert len(x)==30 and len(y)==34
        for xx,yy in zip(x,y[:30]):
            if xx!=yy:
                assert xx.strip()=='' and yy=='NA',f'Value mismatch on row {n+2}'
                changes+=1
        n+=1
Path('work/cod_pair_integrity.json').write_text(json.dumps(dict(rows=n,columns_checked=30,all_values_and_order_match=True,blank_to_NA_changes=changes),indent=2))
print('PASS original 30 columns and record order:',n,'rows; blank-to-NA changes:',changes)
