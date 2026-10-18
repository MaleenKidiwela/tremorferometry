import json, os
v = "nc"; o = "/home/jovyan/tremorferometry/lfe_extension/results/refine"
cfg = dict(variant=v, prod="prod", prod_S13="prod_S13",
           truth={"2010": {"S9": f"r3t2010_{v}"},
                  "2011": {"S15": f"r3t2011_{v}", "S14": f"t2011_S14_{v}", "S13": f"t2011_S13_{v}"},
                  "2012": {"S15": f"r3t2012_{v}", "S14": f"t2012_S14_{v}", "S13": f"t2012_S13_{v}"},
                  "2013": {"S15": f"t2013_S15_{v}", "S14": f"t2013_S14_{v}", "S13": f"t2013_S13_{v}"}})
os.makedirs(f"{o}/conf", exist_ok=True); json.dump(cfg, open(f"{o}/conf/runs.json", "w"), indent=1)
