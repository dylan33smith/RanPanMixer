#!/usr/bin/env bash
# PanMixer's own knob, swept. COMPARISON ARM ONLY.
#
# `capacity` is its utility-loss budget as a fraction 0-1 (obfuscate.py:10). It is
# the dial that decides how much of the target's genotype column gets rewritten, so
# it is the axis to align against our `tau` on a target_fidelity frontier.
#
# ⚠ This runs PanMixer's pinned checkout on PanMixer's own data, unmodified. It is a
# measurement of their tool, never an input to ours; nothing in src/ configs/ tests/
# reads anything it produces.
set -euo pipefail
export PATH=/home/ds85/miniconda3/envs/panmixer/bin:$PATH
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
cd /data/ds85/RanPanMixer/external/PanMixer
OUT=/data/ds85/RanPanMixer/runs/capacity_sweep
mkdir -p "$OUT"
for cap in 0.01 0.05 0.1 0.2 0.3 0.5 0.7 0.9; do
  tag="cap$(echo "$cap" | tr . p)"
  echo "### capacity=$cap  $(date +%H:%M:%S)"
  python tools/panmixer/obfuscate.py starting_data HG00438 "$cap" 21 "$OUT/$tag" --seed 123
done
echo "### SWEEP DONE $(date +%H:%M:%S)"
