import os
import csv, gzip, pathlib, collections

root = pathlib.Path(__file__).resolve().parent.parent
src = pathlib.Path(os.environ["SEER_COD_EXPORT"])

cardio = {
    "Ischemic heart disease", "Cerebrovascular diseases",
    "Other and unspecified disorders of the circulatory system",
    "Hypertensive disease", "Diseases of arteries, arterioles and capillaries",
    "Pulmonary heart disease and diseases of pulmonary circulation",
}
respiratory = {
    "Pneumonia and Influenza",
    "Chronic Obstructive Pulmonary Disease and Allied Conditions",
}
unknown = {"State DC not available or state DC available but no COD", "NA", ""}
noncancer = cardio | respiratory | unknown | {
    "Alzheimers (ICD-9 and ICD-10 only)",
    "Other infectious and Parasitic Diseases incl HIV", "Other COD",
    "Chronic Liver Disease and Cirrhosis", "Septicemia",
    "Stomach and Duodenal Ulcers", "Accidents and Adverse Effects",
    "Nephritis, Nephrotic Syndrome and Nephrosis", "COVID (2020+ only)",
    "Congenital Anomalies", "Diabetes Mellitus",
    "Symptoms, Signs and Ill-Defined Conditions", "Suicide and Self-Inflicted Injury",
    "Complications of medical and surgical care (Y40-Y84, Y88) (ICD-10 only, 1999+)",
    "Homicide and Legal Intervention", "Complications of Pregnancy, Childbirth, Puerperium",
}
nonmalignant_neoplasm = {
    "Benign and Borderline: All Other sites",
    "Brain, CNS Other and Intracranial Gland (Benign and Borderline)",
    "In situ neoplasms",
}

def detail_category(event, cod):
    if event == 0: return "alive"
    if event == 1: return "breast"
    if event == 3 or cod in unknown: return "unknown"
    if cod in cardio: return "cardiovascular"
    if cod in respiratory: return "respiratory"
    if cod == "COVID (2020+ only)": return "covid"
    if cod not in noncancer and cod not in nonmalignant_neoplasm and cod != "Breast":
        return "other_cancer"
    return "remaining_other"

counts = collections.Counter()
fields = ["row","sex","year","age","race","stage","grade","er","months","days",
          "event","detail","cod","followup_year","hist","mainhist"]
with src.open(encoding="utf-8-sig", newline="") as f, gzip.open(root/"work/cod_extended.csv.gz", "wt", encoding="utf-8", newline="") as o:
    rd = csv.DictReader(f); wr = csv.DictWriter(o, fieldnames=fields); wr.writeheader()
    for line, r in enumerate(rd, 2):
        if r["Type of Reporting Source"] in {"Death certificate only", "Autopsy only"}: continue
        try: months = int(r["Survival months"])
        except ValueError: continue
        if months >= 9999: continue
        if r["Sequence number"] not in {"One primary only", "1st of 2 or more primaries"}: continue
        year = int(r["Year of diagnosis"]); hist = int(r["Histologic Type ICD-O-3"])
        st = r["Summary stage 2000 (1998-2017)"] if year < 2004 else r["Combined Summary Stage with Expanded Regional Codes (2004+)"]
        stage = next((x for x in ("Localized","Regional","Distant") if st.startswith(x)), "Unknown")
        raw = r["Grade Recode (thru 2017)"] if year < 2018 else r["Derived Summary Grade 2018 (2018+)"]
        gm = {"Well differentiated; Grade I":"I","Moderately differentiated; Grade II":"II","Poorly differentiated; Grade III":"III","Undifferentiated; anaplastic; Grade IV":"III","Site-specific grade system category (1)":"I","Site-specific grade system category (2)":"II","Site-specific grade system category (3)":"III","Well differentiated":"I","Moderately differentiated":"II","Poorly differentiated":"III","Undifferentiated and anaplastic":"III"}
        grade = gm.get(raw, "Unknown")
        er = r["ER Status Recode Breast Cancer (1990+)"]
        if er not in {"Positive","Negative"}: er = "Unknown"
        race = {"Non-Hispanic White":"NHW","Non-Hispanic Black":"NHB","Hispanic (All Races)":"Hispanic","Non-Hispanic Asian or Pacific Islander":"NHAPI","Non-Hispanic American Indian/Alaska Native":"Other"}.get(r["Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)"], "Other")
        dead = r["Vital status recode (study cutoff used)"] == "Dead"
        bc = r["SEER cause-specific death classification"] == "Dead (attributable to this cancer dx)"
        other = r["SEER other cause of death classification"] == "Dead (attributable to causes other than this cancer dx)"
        event = 0 if not dead else 1 if bc else 2 if other else 3
        cod = r["COD to site recode ICD-O-3 2023 Revision Expanded (1999+)"]
        detail = detail_category(event, cod)
        rec = dict(row=line, sex=r["Sex"], year=year, age=int(r["Age recode with single ages and 85+"].split()[0].replace("+","")), race=race, stage=stage, grade=grade, er=er, months=months, days=r["Survival Days"], event=event, detail=detail, cod=cod, followup_year=r["Year of follow-up recode"], hist=hist, mainhist=int(8010 <= hist <= 8589))
        wr.writerow(rec); counts[(r["Sex"], detail, year <= 2018, rec["mainhist"])] += 1
print(counts)
