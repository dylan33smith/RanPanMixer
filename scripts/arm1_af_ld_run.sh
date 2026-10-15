#!/bin/bash
# One PanMixer in-graph utility evaluator over all 308 rows. Completeness asserted from disk.
set -u
R=/data/ds85/RanPanMixer/runs/beagle_attack/root
L=/data/ds85/RanPanMixer/runs/beagle_attack/logs
PY=/home/ds85/miniconda3/envs/panmixer/bin/python
NAME=$1; NPAR=$2
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
run_one() {
  local row=$1 cwd="/data/ds85/RanPanMixer/runs/beagle_attack/scratch/${NAME}_row${row}"
  mkdir -p "$cwd"
  ( cd "$cwd" && PYTHONPATH="$R" PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg \
      "$PY" -B "$R/tools/downstream/utility_in/${NAME}.py" 5 False 21 "$row" \
      > "$L/${NAME}_row${row}.log" 2>&1 )
  echo $? > "$L/${NAME}_row${row}.status"
}
for row in $(seq 0 307); do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 2; done
  run_one "$row" &
done
wait
ok=0; bad=0
for row in $(seq 0 307); do
  [ "$(cat "$L/${NAME}_row${row}.status" 2>/dev/null)" = "0" ] && ok=$((ok+1)) || { bad=$((bad+1)); echo "FAILED row $row"; }
done
echo "$NAME: 308 launched, $ok ok, $bad failed"
[ "$bad" -eq 0 ]
