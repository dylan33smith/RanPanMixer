#!/bin/bash
# One PanMixer in-graph utility evaluator over every manifest row. Completeness asserted
# from disk. ROOT/LOGS/SCRATCH are overridable so a second mirror never touches this one.
#
# ⚠⚠ EXP AND ROWSKEY ARE IN THE LOG AND SCRATCH NAMES ON PURPOSE. The 2026-10-28 our-build
# run evaluates TWO experiments from one root -- exp_5 (ours, `rows`) and exp_6 (theirs,
# `theirs_rows`) -- and both start at row 0. Without EXP in the name, `af_loss_row0.status`
# would be written twice and the completeness count would be over one experiment while
# claiming both. Defaults reproduce the pre-2026-10-28 single-experiment behaviour.
set -u
R=${ROOT:-/data/ds85/RanPanMixer/runs/beagle_attack/root}
L=${LOGS:-/data/ds85/RanPanMixer/runs/beagle_attack/logs}
SCR=${SCRATCH:-/data/ds85/RanPanMixer/runs/beagle_attack/scratch}
EXP=${EXP:-5}
ROWSKEY=${ROWSKEY:-rows}
PY=/home/ds85/miniconda3/envs/panmixer/bin/python
NAME=$1; NPAR=$2
N=$(PATH=/home/ds85/miniconda3/envs/ranpanmixer/bin:$PATH python -c \
   "import json;print(len(json.load(open('$R/manifest.json'))['$ROWSKEY']))")
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
TAG="${NAME}_exp${EXP}"
run_one() {
  local row=$1 cwd="$SCR/${TAG}_row${row}"
  mkdir -p "$cwd"
  ( cd "$cwd" && PYTHONPATH="$R" PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg \
      "$PY" -B "$R/tools/downstream/utility_in/${NAME}.py" "$EXP" False 21 "$row" \
      > "$L/${TAG}_row${row}.log" 2>&1 )
  echo $? > "$L/${TAG}_row${row}.status"
}
for row in $(seq 0 $((N-1))); do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 2; done
  run_one "$row" &
done
wait
ok=0; bad=0
for row in $(seq 0 $((N-1))); do
  [ "$(cat "$L/${TAG}_row${row}.status" 2>/dev/null)" = "0" ] && ok=$((ok+1)) || { bad=$((bad+1)); echo "FAILED row $row"; }
done
echo "$TAG: $N launched, $ok ok, $bad failed"
[ "$bad" -eq 0 ]
