#!/usr/bin/env bash
# P3 templates: (1) freeze the G3 era templates (POLARIS + LZB, 2004 build blocks) for production; (2) SNB-CMG40T A/B:
# SNB era template (2004 build blocks, re-timed without SNB) vs the 2012 SNB template inside the era detector, scored on the
# already-spent 2004-07 TEST blocks with the G3 scorer.  Rule (fixed here, before running): adopt the SNB era template for
# SNB days before 2006-04-20 iff its held-out recall is >= the 2012-template recall at both matched scrambled-control rates
# AND the rediscovered-family count does not fall.  Writes templates/p3_config.json and P3_DONE / P3_FAILED.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
echo $$ > $O/P3_PID
GO=/home/jovyan/tremorferometry/svi_science/results/backward_gates
fail(){ echo "FAILED $1 $(date -u +%FT%TZ)"; date -u +%FT%TZ > $O/P3_FAILED; exit 1; }
cp -n $GO/tables/g3_era_templates.npz $GO/tables/g3_era_weights.csv $GO/tables/g3_era_ampref.csv $GO/tables/g3_carriers.csv $O/templates/ || fail copy
sha256sum $O/templates/g3_era_* > $O/templates/SHA256SUMS
OMP_NUM_THREADS=4 python $S/bw_prod_p3_snb.py > $O/logs/p3_snb.log 2>&1 || fail snb_build
export SVI_REF_DIR=$GO/mf BW_SITELIST_EXT=1
for i in 0 1; do
  BW_ERA_TPL=$O/templates/g3_era_templates.npz BW_TPL_OVERRIDE=$O/templates/snb40t_2004.npz \
    python $S/bw_gate_mf.py --backward -- --run p3c_CN4ERA_SNB --kind truth --year 2012 --sset CN4ERA --clean 0 \
    --days-file $GO/tables/g3_days.json --cache-dirs $CB --shard $i --nshard 2 > $O/logs/p3c_s$i.log 2>&1 &
done
wait
BW_G3_RUNS=g3b_CN4ERA,p3c_CN4ERA_SNB BW_G3_TAG=p3_snb python $S/bw_gate_g3_score.py > $O/logs/p3_score.log 2>&1 || fail score
python - <<PY
import json, shutil
d = json.load(open("$GO/tables/p3_snb_decision.json"))
rows = d["rows"]
adopt = all(r["recall_b"] >= r["recall_a"] and r["fam_red_b"] >= r["fam_red_a"] for r in rows)
cfg = dict(era_npz="$O/templates/g3_era_templates.npz", snb_override="$O/templates/snb40t_2004.npz" if adopt else "",
           snb_override_until="2006-04-20", snb_ab=[{k: r[k] for k in ("T_ref", "recall_a", "recall_b", "fam_red_a", "fam_red_b")} for r in rows],
           rule="adopt SNB era template iff recall >= 2012-template recall at both matched rates and families not fewer (2004-07 test blocks)")
json.dump(cfg, open("$O/templates/p3_config.json", "w"), indent=1)
for f in ("p3_snb_recall.csv", "p3_snb_decision.json", "p3_snb_family_recall.csv", "p3_snb_quiet.csv"):
    shutil.copy("$GO/tables/" + f, "$O/tables/" + f)
print("SNB era template adopted:", adopt)
PY
date -u +%FT%TZ > $O/P3_DONE
