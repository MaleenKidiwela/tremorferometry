#!/bin/bash
# Track A step 0(a)+(b) -- re-score every E6 pass with --e6 (fixes Merlin F4: the E5r key L{L}_T{int(T)} collides for the
# half-integer thresholds, so the on-disk e6_score_e6_{sites,mergeA,mergeB}.json label 7.5 as "T7" and 8.5 as "T8").
# Inputs stay in lfe_catalog_repair/results (read-only); outputs go to lfe_extension/results/e6_rescore/.
# $1 = "" (nearest-neighbour, reproduces E6) or "--one2one" (F11 greedy-by-time one-to-one); $2 = output sub-directory.
set -u
cd /home/jovyan/tremorferometry/lfe_extension
IN=/home/jovyan/tremorferometry/lfe_catalog_repair/results
M=${1:-}; SUB=${2:-nn}
OUT=results/e6_rescore/$SUB; mkdir -p $OUT logs
run () { tag=$1; shift
  python scripts/e5r_score.py --tag $tag --e6 --in-dir $IN --out-dir $OUT $M "$@" > logs/a0_score_${SUB}_${tag}.log 2>&1 \
    && echo "  OK   $tag" || { echo "  FAIL $tag"; tail -5 logs/a0_score_${SUB}_${tag}.log; }
}
ST15=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
ST13=PGC,NLLB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
ST8=PGC,VGZ,NLLB,SNB,B009,B010,B011,B926
echo "[$(date -u +%FT%TZ)] rescore ($SUB, matching='${M:-nearest}')"
run e6_sites          --qc-csv $IN/e5r_template_qc_2012.csv   --stations $ST15
run e6_mergeA         --qc-csv $IN/e5r_template_qc_2012.csv   --stations $ST15
run e6_mergeB         --qc-csv $IN/e5r_template_qc_2012.csv   --stations $ST15
run e6_mergeB_cap25   --qc-csv $IN/e5r_template_qc_2012.csv   --stations $ST15
run e6_mergeB_coinc   --qc-csv $IN/e5r_template_qc_2012.csv   --stations $ST15
run e6_lev_a          --qc-csv results/e6_rescore/qc/qc_from_npz_lev_a.csv  --stations $ST15
run e6_lev_b          --qc-csv $IN/e5r_template_qc_2012.csv   --stations $ST8
run e6_lev_c          --qc-csv results/e6_rescore/qc/qc_from_npz_lev_c.csv  --stations $ST15
run e6_1b_noVGZSNB    --qc-csv $IN/e5r_template_qc_2012.csv   --stations $ST13
echo "[$(date -u +%FT%TZ)] done -> $OUT"
