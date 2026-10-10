#!/usr/bin/env bash
# Build the cohort artifacts for all 22 autosomes. The driver docs/plan.md has been
# calling for since 2026-09-29.
#
# ⚠ FAN-OUT IS BOUNDED BY MEMORY, NOT BY CORES. vg deconstruct peaked at 9.48 GB RSS
# on chr21's 435 MB GFA -- about 22 MB of RSS per MB of gzipped graph -- so chr1's
# 1.2 GB GFA implies roughly 26 GB. Running 22 at once would ask for several hundred
# GB against 251. The jobs are therefore ordered LARGEST FIRST and run at a fan-out
# computed from the actual file sizes, so the biggest chromosomes do not land together.
#
# ⚠ AND IT ASSERTS ITS OWN COMPLETENESS. A driver that filters results cannot tell
# "found nothing" from "ran nothing" (CLAUDE.md), and a partial genome looks exactly
# like a whole one to everything downstream.
set -uo pipefail
cd "$(dirname "$0")/.."

FANOUT="${FANOUT:-3}"
LOGDIR="logs/cohorts"; mkdir -p "$LOGDIR"

# largest GFA first, so the memory peaks are spread rather than clustered
mapfile -t ORDER < <(for c in $(seq 1 22); do
  f="data/raw/chr${c}.hprc-v1.0-pggb.gfa.gz"
  [ -s "$f" ] && echo "$(stat -c %s "$f") $c"
done | sort -rn | awk '{print $2}')

echo "order (largest GFA first): ${ORDER[*]}"
echo "fan-out: $FANOUT"

run_one () {
  local c="$1"
  local cfg="configs/inputs.chr${c}.yaml"
  if [ ! -f "$cfg" ]; then echo "chr$c: NO CONFIG" > "$LOGDIR/chr${c}.log"; return 1; fi
  prep-cohort --chrom "chr${c}" --config "$cfg" > "$LOGDIR/chr${c}.log" 2>&1
  local rc=$?
  echo "$rc" > "$LOGDIR/chr${c}.status"
  return $rc
}
export -f run_one; export LOGDIR

printf '%s\n' "${ORDER[@]}" | xargs -P "$FANOUT" -I{} bash -c 'run_one {}'

# --- completeness, asserted ---
echo; echo "=== completeness ==="
bad=0
for c in $(seq 1 22); do
  st="$LOGDIR/chr${c}.status"
  d="data/cohort/chr${c}"
  if [ ! -s "$st" ]; then echo "chr$c: NEVER RAN"; bad=$((bad+1)); continue; fi
  if [ "$(cat "$st")" != "0" ]; then echo "chr$c: EXIT $(cat "$st")"; bad=$((bad+1)); continue; fi
  for f in haplotypes.npy positions.npy support.npy genetic_pos.npy in_chain.npy \
           cm_source.npy sites.tsv manifest.json; do
    [ -s "$d/$f" ] || { echo "chr$c: missing $f"; bad=$((bad+1)); }
  done
done
if [ "$bad" -ne 0 ]; then
  echo; echo "INCOMPLETE: $bad problem(s). Do NOT quote a genome-wide number from this tree."
  exit 1
fi
echo "all 22 autosomes built and complete"
