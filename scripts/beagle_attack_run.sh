#!/bin/bash
# PanMixer's BEAGLE RECONSTRUCTION ATTACK on all 308 of our ARM 1 releases.
#
# ⚠ Every scoring step is THEIR code, called not copied (see scripts/beagle_attack_one.py).
#   This file only schedules it.
# ⚠ Nothing under external/ is written: the mirror root holds a COPY of constants.py with
#   tools/, starting_data/ and downloaded_tools/ symlinked back read-only.
# ⚠ PYTHONDONTWRITEBYTECODE=1 is LOAD-BEARING: their __pycache__ holds a 3.12 .pyc and the
#   panmixer env is 3.11, so an import would otherwise write a new .pyc under theirs.
# ⚠ Each row gets its own CWD because their code writes some outputs relative to the working
#   directory; concurrent rows from one CWD would race.
# ⚠ Completeness is asserted FROM DISK at the end, never from this loop's own bookkeeping.
#
#   NPAR=8 THREADS=6 bash scripts/beagle_attack_run.sh            # their protocol
#   MODE=selfout ROWS="0 6 7 13" bash scripts/beagle_attack_run.sh # the labelled diagnostic
set -u
REPO=/home/ds85/projects/Pangenome/RanPanMixer
ROOT=/data/ds85/RanPanMixer/runs/beagle_attack/root
LOGS=/data/ds85/RanPanMixer/runs/beagle_attack/logs
PY=/home/ds85/miniconda3/envs/panmixer/bin/python
NPAR=${NPAR:-8}
THREADS=${THREADS:-6}
XMX=${XMX:-12g}
MODE=${MODE:-protocol}
N=$(PATH=/home/ds85/miniconda3/envs/ranpanmixer/bin:$PATH python -c \
   "import json;print(len(json.load(open('$ROOT/manifest.json'))['rows']))")
ROWS=${ROWS:-$(seq 0 $((N-1)))}
mkdir -p "$LOGS"

EXTRA=""
SUF=""
if [ "$MODE" = "selfout" ]; then EXTRA="--exclude-self"; SUF="_selfout"; fi

run_one() {
  local row=$1
  local tag="beagle${SUF}_row${row}"
  local cwd="/data/ds85/RanPanMixer/runs/beagle_attack/scratch/$tag"
  mkdir -p "$cwd"
  ( cd "$cwd" \
    && PYTHONDONTWRITEBYTECODE=1 "$PY" -B "$REPO/scripts/beagle_attack_one.py" \
         --row "$row" --threads "$THREADS" --xmx "$XMX" $EXTRA \
       > "$LOGS/$tag.log" 2>&1 )
  echo $? > "$LOGS/$tag.status"
}

echo "BEAGLE RECONSTRUCTION ATTACK on our arm 1 releases"
echo "  mode: $MODE   rows: $(echo $ROWS | wc -w)   parallel: $NPAR x $THREADS threads, -Xmx$XMX"
echo "  44 targets x 7 tau, chr21, K = 86 release scored on their K = 88 axis"
echo "  ⚠ 39 of 44 targets are IN the reference panel; their protocol excludes nobody"
echo
started=0
for row in $ROWS; do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 5; done
  run_one "$row" &
  started=$((started+1))
  if [ $((started % 25)) -eq 0 ]; then echo "  launched $started..."; fi
done
wait

ok=0; bad=0
for row in $ROWS; do
  f="$LOGS/beagle${SUF}_row${row}.status"
  if [ -f "$f" ] && [ "$(cat "$f")" = "0" ]; then ok=$((ok+1));
  else bad=$((bad+1)); echo "  FAILED: row $row -> $(cat "$f" 2>/dev/null || echo no-status)"; fi
done
echo
echo "$started launched, $ok succeeded, $bad failed"
RESNAME="accuracy_results.json"
[ "$MODE" = "selfout" ] && RESNAME="accuracy_results_selfout.json"
present=$(find "$ROOT/experiments" -name "$RESNAME" | wc -l)
echo "result files on disk: $present (expect $(echo $ROWS | wc -w))"
[ "$bad" -eq 0 ]
