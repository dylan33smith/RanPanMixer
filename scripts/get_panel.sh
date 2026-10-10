#!/usr/bin/env bash
# Fetch the 1000G 30x GRCh38 phased panel for chr1-chr20 and chr22.
#
# ⚠ WHY THIS PANEL AND NOT PHASE 3. The pangenome is GRCh38 and this panel is
# GRCh38-called. PanMixer's released pipeline wires in 1000G Phase 3, which is
# GRCh37, and the mismatch is measurable: 0.30% of chr21 pangenome records match
# Phase 3 on (POS, REF, ALT) against 75.88% here -- a 252x difference attributable to
# the build alone. Decision F is our commitment not to inherit that.
#
# It serves two purposes and they must not be conflated: it is where Arm 2's external
# TARGETS come from, and it is the ATTACK DATABASE. docs/data.md carries the
# eligibility rules that separate the two.
#
# ⚠ Writes a sha256 beside each file. An artifact with no recorded checksum is not
# quotable, and re-deriving one later means re-reading 28 GB.
set -uo pipefail
B="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
STEM="1kGP_high_coverage_Illumina"
SUF="filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
OUT="${1:?usage: get_panel.sh <outdir>}"
mkdir -p "$OUT"
# ⚠ PARALLEL, because a single connection to EBI measured about 700 KB/s, which
# would be 11 hours for 28 GB. Four at once is polite and roughly four times faster;
# the limit is per-connection throughput, not our bandwidth.
# ⚠ And `curl -C -` RESUMES a partial file, so a restart or a dropped connection does
# not re-download what is already on disk.
export B STEM SUF OUT
fetch_one () {
  local c="$1" ext="$2"
  local f="$OUT/1000g_30x_chr${c}.vcf.gz${ext}"
  [ -s "$f" ] && return 0
  local url="$B/${STEM}.chr${c}.${SUF}${ext}"
  if curl -sSfL -C - --retry 5 --retry-delay 5 -o "$f.part" "$url"; then
    mv "$f.part" "$f"
    [ -z "$ext" ] && sha256sum "$f" | awk '{print $1}' > "$f.sha256"
    echo "ok chr$c$ext  $(stat -c %s "$f") B"
  else
    echo "FAILED chr$c$ext"; return 1
  fi
}
export -f fetch_one
jobs_list=""
for c in $(seq 1 22); do
  [ "$c" = "21" ] && continue          # already on disk as 1000g_30x_phased.vcf.gz
  jobs_list="$jobs_list $c:  $c:.tbi"
done
printf '%s\n' $jobs_list | tr ' ' '\n' | grep -v '^$' \
  | xargs -P 4 -I{} bash -c 'IFS=: read -r c e <<< "{}"; fetch_one "$c" "$e"'
fail=0
# ⚠ Assert completeness and fail loudly: a partial panel is a partial target pool
# AND a partial attack database, and neither announces itself.
missing=0
for c in $(seq 1 22); do
  [ "$c" = "21" ] && continue
  for ext in "" ".tbi"; do
    [ -s "$OUT/1000g_30x_chr${c}.vcf.gz${ext}" ] || { echo "MISSING chr$c$ext"; missing=$((missing+1)); }
  done
done
echo "failures=$fail missing=$missing"
[ "$fail" -eq 0 ] && [ "$missing" -eq 0 ]
