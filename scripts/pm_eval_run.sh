#!/bin/bash
# Run PanMixer's privacy evaluators, UNMODIFIED, on both tools' releases.
#
# ⚠⚠ Nothing under external/ is written. The mirror root holds a COPY of constants.py,
#    so BASE_PATH resolves there, while tools/ and starting_data/ are read-only symlinks.
# ⚠  PYTHONDONTWRITEBYTECODE=1 is LOAD-BEARING: their __pycache__ holds a 3.12 .pyc and
#    the panmixer env is 3.11, so an import would otherwise write a new .pyc under theirs.
# ⚠  Each run gets its OWN CWD because gap_score.py:74 and diploid_gap_score.py:100 call
#    plt.savefig("test_hist.png") relative to the working directory; concurrent runs from
#    one CWD would race on that file. BASE_PATH is unaffected -- constants.py derives it
#    from __file__, not from CWD.
# ⚠  One experiment per evaluator: gap_score and diploid_gap_score write IDENTICAL output
#    filenames into the row directory, so sharing a row would let one overwrite the other.
set -u
ROOT=/data/ds85/RanPanMixer/runs/pm_eval
PY=/home/ds85/miniconda3/envs/panmixer/bin/python
LOGS=$ROOT/logs
NPAR=${NPAR:-6}
mkdir -p "$LOGS"

run_one() {
  local exp=$1 name=$2 row=$3
  local tag="${name}_row${row}"
  local cwd="$ROOT/scratch/$tag"
  mkdir -p "$cwd"
  ( cd "$cwd" \
    && PYTHONPATH="$ROOT" PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg \
       "$PY" -B "$ROOT/tools/downstream/privacy/${name}.py" "$exp" 21 "$row" \
       > "$LOGS/$tag.log" 2>&1 )
  echo $? > "$LOGS/$tag.status"
}

echo "PanMixer privacy evaluators on BOTH tools' releases"
echo "  rows 0-6   = OURS   (tau 0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)"
echo "  rows 7-14  = THEIRS (capacity 0.01 .. 0.9)"
echo "  evaluators: gap_score (exp 1), diploid_gap_score (exp 2), MIA_privacy (exp 3)"
echo "  ⚠ sampled_gap_score is NOT runnable: its pre-stacker input"
echo "    starting_data/chr21/subjects/ does not exist."
echo "  parallelism: $NPAR"
echo

started=0
for spec in "1 gap_score" "2 diploid_gap_score" "3 MIA_privacy"; do
  set -- $spec; exp=$1; name=$2
  for row in $(seq 0 14); do
    while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 3; done
    run_one "$exp" "$name" "$row" &
    started=$((started+1))
  done
done
wait

ok=0; bad=0
for f in "$LOGS"/*.status; do
  if [ "$(cat "$f")" = "0" ]; then ok=$((ok+1)); else bad=$((bad+1)); echo "  FAILED: $(basename "$f" .status) -> $(cat "$f")"; fi
done
echo
echo "$started launched, $ok succeeded, $bad failed"
# ⚠ assert completeness from the DISK, not from the loop's own bookkeeping.
present=$(find "$ROOT/experiments" -name "gap_score_all.json" -o -name "MIA_privacy.json" | wc -l)
echo "result files on disk: $present (expect 45: 15 rows x 3 evaluators)"
[ "$bad" -eq 0 ] || exit 1
