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
# ⚠⚠ THREE WAYS THIS USED TO REPORT SUCCESS HAVING RUN NOTHING, all found by audit
# 2026-10-28. (1) On a ROWSKEY that is not a manifest key the python one-liner raised,
# the command substitution yielded "", the ASSIGNMENT still succeeded so `set -u` never
# fired, `seq 0 $((N-1))` became `seq 0 -1` = no jobs, and the tally printed
# "0 launched, 0 ok, 0 failed" and exited 0 -- which on a re-run tables the PREVIOUS
# run's outputs as current. (2) EXP and ROWSKEY were independent, so EXP=6 ROWSKEY=rows
# counted 308 of exp_6's 352 rows and called it complete, leaving ~5.5 individuals out
# of every column under an "n = 44" header. (3) A bad key is now fatal here rather than
# arithmetic later.
N=$(PATH=/home/ds85/miniconda3/envs/ranpanmixer/bin:$PATH python -c \
   "import json,sys
m=json.load(open('$R/manifest.json'))
k='$ROWSKEY'
if k not in m: sys.exit(f'manifest has no key {k!r}; keys: {sorted(m)}')
rows=m[k]
exp='$EXP'
want={'rows':'5','theirs_rows':'6'}
if want.get(k) not in (None, exp): sys.exit(f'ROWSKEY {k!r} belongs to exp_{want[k]}, not exp_{exp}')
if not rows: sys.exit(f'{k} is empty')
print(len(rows))") || exit 1
case "$N" in ''|*[!0-9]*) echo "FATAL: row count not a positive integer: '$N'" >&2; exit 1;; esac
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
