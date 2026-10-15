#!/bin/bash
# PanMixer's three runnable LINKAGE/MIA evaluators on all 308 of our ARM 1 releases.
#
# ⚠ Their code, UNMODIFIED. Nothing under external/ is written: the mirror root holds a COPY
#   of constants.py with tools/ and starting_data/ symlinked back read-only.
# ⚠ PYTHONDONTWRITEBYTECODE=1 is LOAD-BEARING: their __pycache__ holds a 3.12 .pyc and the
#   panmixer env is 3.11, so an import would otherwise write a new .pyc under theirs.
# ⚠ Each run gets its OWN CWD because gap_score.py:74 and diploid_gap_score.py:100 call
#   plt.savefig("test_hist.png") relative to the working directory; concurrent runs from one
#   CWD would race on that file. BASE_PATH is unaffected -- constants.py derives it from
#   __file__, not from CWD.
# ⚠ MEMORY IS THE BINDING CONSTRAINT, not cores: each job loads 1000g_phased_masked.npy
#   (3.1 GB) and diploid_gap_score also loads 1000g_phased_masked_posref.npy (3.2 GB). NPAR
#   is therefore low by default and must stay that way.
# ⚠ Completeness is asserted FROM DISK at the end, never from this loop's own bookkeeping.
#
#   NPAR=8 bash scripts/arm1_pm_attacks_run.sh
set -u
ROOT=/data/ds85/RanPanMixer/runs/beagle_attack/root
LOGS=/data/ds85/RanPanMixer/runs/beagle_attack/logs
PY=/home/ds85/miniconda3/envs/panmixer/bin/python
NPAR=${NPAR:-8}
N=$(PATH=/home/ds85/miniconda3/envs/ranpanmixer/bin:$PATH python -c \
   "import json;print(len(json.load(open('$ROOT/manifest.json'))['rows']))")
ROWS=${ROWS:-$(seq 0 $((N-1)))}
mkdir -p "$LOGS"

run_one() {
  local exp=$1 name=$2 row=$3
  local tag="${name}_row${row}"
  local cwd="/data/ds85/RanPanMixer/runs/beagle_attack/scratch/$tag"
  mkdir -p "$cwd"
  ( cd "$cwd" \
    && PYTHONPATH="$ROOT" PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg \
       "$PY" -B "$ROOT/tools/downstream/privacy/${name}.py" "$exp" 21 "$row" \
       > "$LOGS/$tag.log" 2>&1 )
  echo $? > "$LOGS/$tag.status"
}

echo "PanMixer's linkage/MIA evaluators on our ARM 1 releases"
echo "  308 rows = 44 targets x 7 tau; 3 evaluators = $((N*3)) jobs"
echo "  exp_2 gap_score, exp_3 diploid_gap_score, exp_4 MIA_privacy"
echo "  ⚠ sampled_gap_score is NOT runnable: its pre-stacker input does not exist here"
echo "  parallelism: $NPAR"
echo
started=0
for spec in "2 gap_score" "3 diploid_gap_score" "4 MIA_privacy"; do
  set -- $spec; exp=$1; name=$2
  echo "--- $name (exp_$exp) ---"
  for row in $ROWS; do
    while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 4; done
    run_one "$exp" "$name" "$row" &
    started=$((started+1))
    if [ $((started % 50)) -eq 0 ]; then echo "  launched $started..."; fi
  done
  wait
done

ok=0; bad=0
for spec in "2 gap_score" "3 diploid_gap_score" "4 MIA_privacy"; do
  set -- $spec; name=$2
  for row in $ROWS; do
    f="$LOGS/${name}_row${row}.status"
    if [ -f "$f" ] && [ "$(cat "$f")" = "0" ]; then ok=$((ok+1));
    else bad=$((bad+1)); echo "  FAILED: ${name} row $row -> $(cat "$f" 2>/dev/null || echo no-status)"; fi
  done
done
echo
echo "$started launched, $ok succeeded, $bad failed"
present=$(find "$ROOT/experiments/exp_2" "$ROOT/experiments/exp_3" -name "gap_score_all.json" | wc -l)
presentm=$(find "$ROOT/experiments/exp_4" -name "MIA_privacy.json" | wc -l)
echo "result files on disk: gap $present (expect $((N*2))), MIA $presentm (expect $N)"
[ "$bad" -eq 0 ]
