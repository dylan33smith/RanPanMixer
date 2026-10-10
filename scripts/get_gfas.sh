#!/usr/bin/env bash
# Fetch the 21 remaining HPRC v1.0 PGGB chromosome GFAs (decision D: autosomes only).
# ⚠ Writes a sha256 beside each file. An artifact with no recorded checksum is not
# quotable under docs/data.md's rules, and re-deriving it later means re-reading 14 GB.
set -uo pipefail
B="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms"
OUT="${1:?usage: get_gfas.sh <outdir>}"
mkdir -p "$OUT"
fail=0
for c in $(seq 1 22); do
  [ "$c" = "21" ] && continue          # already on disk
  f="$OUT/chr${c}.hprc-v1.0-pggb.gfa.gz"
  if [ -s "$f" ] && [ -s "$f.sha256" ]; then echo "skip chr$c (have it)"; continue; fi
  echo "--> chr$c"
  if curl -sSfL --retry 3 --retry-delay 5 -o "$f.part" "$B/chr${c}.hprc-v1.0-pggb.gfa.gz"; then
    mv "$f.part" "$f"
    sha256sum "$f" | awk '{print $1}' > "$f.sha256"
    echo "    $(stat -c %s "$f") B  $(cat "$f.sha256")"
  else
    echo "    FAILED chr$c"; rm -f "$f.part"; fail=$((fail+1))
  fi
done
# ⚠ Assert completeness and fail loudly. A loop that filters cannot tell
# "found nothing" from "ran nothing" (CLAUDE.md).
missing=0
for c in $(seq 1 22); do
  [ -s "$OUT/chr${c}.hprc-v1.0-pggb.gfa.gz" ] || { echo "MISSING chr$c"; missing=$((missing+1)); }
done
echo "failures=$fail missing_of_22=$missing"
[ "$fail" -eq 0 ] && [ "$missing" -eq 0 ]
