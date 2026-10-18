import pandas as pd, glob, sys
sta=sys.argv[1:]
rows=[]
for f in glob.glob("mus_*.csv"):
    _,s,c,m=f[:-4].split("_",3)
    if s not in sta: continue
    try: D=pd.read_csv(f,skiprows=1)
    except Exception: continue
    if D.empty: continue
    D["month"]=D.start.str[:7].str.replace("/","-")
    g=D.groupby("month").value.median() if m not in("num_gaps",) else D.groupby("month").value.sum()
    for k,v in g.items(): rows.append((s,c,m,k,v))
R=pd.DataFrame(rows,columns=["sta","cha","m","month","v"])
P=R.pivot_table(index="month",columns=["m","sta","cha"],values="v")
pd.set_option("display.width",250); pd.set_option("display.max_columns",40); pd.set_option("display.max_rows",100)
for m in ["sample_rms","percent_availability","num_gaps","timing_quality","clock_locked","sample_mean"]:
    if m in P.columns.get_level_values(0): print("==",m); print(P[m].round(1).to_string())
