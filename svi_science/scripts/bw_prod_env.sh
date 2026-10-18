# shared settings for the SVI backward production queues (sourced)
S=/home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/backward
LFE=/home/jovyan/tremorferometry/lfe_extension
CB=$LFE/cache20hz_backward
FWD=$LFE/cache20hz_r1,$LFE/cache20hz,/home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz,/home/jovyan/tremorferometry/lfe_picker/cache20hz
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export BW_LEDGERS=$O/tables/ledger_A.csv,$O/tables/ledger_B.csv,/home/jovyan/tremorferometry/svi_science/results/backward_gates/tables/download_ledger.csv
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
CN4=PGC,VGZ,NLLB,SNB; PO=TWKB,TSJB,KLNB,SILB,SSIB,TWGB; C8=SHVB,MGCB,TWBB,SHDB,KHVB
PB=B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
# retry a resumable python step up to 5 times (download restarts skip finished days)
retry(){ local k; for k in 1 2 3 4 5; do "$@" && return 0; echo "[retry $k] rc=$? $(date -u +%FT%TZ): $*"; sleep 60; done; return 1; }
# --- 2026-10-09 (user): KELB added; NRCan queue reordered: Bostock ETS windows +-2 months first (truth_episodes.csv)
ERA_TPL=$O/templates/kelb/g3_era_templates.npz
# name start end (2003 and 2010-2013 windows are already complete)
WINDOWS="W2004:2004-05-02:2004-09-27 W2005:2005-07-03:2005-11-25 W2007:2006-11-24:2007-04-03 W2008:2008-03-10:2008-08-02 W2009:2009-03-10:2009-08-07"
