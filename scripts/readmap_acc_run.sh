#!/bin/bash
# READ MAPPING ON THE RELEASED GRAPH, on the REPORTED `acc_ceiling` 0.50-0.80 GRID.
#
# EIGHT graphs, ONE code path, so the obfuscated numbers and the baseline are
# apples-to-apples: an unobfuscated BASELINE plus HG00438's release at each of the seven
# ceiling points. The published absolute level is NOT the comparison (their vg is 1.36,
# ours 1.68, and that offset was measured at +1.2/+2.2 pts in September and is
# unexplained after the build hypothesis was refuted). The comparison is OUR delta
# against THEIR delta, which is ~0.01 pt.
#
# ⚠⚠ PRODUCTIONISED FROM A SCRATCHPAD. The 2026-10-26 eta-grid run was driven by
# `readmap_eta.sh`, which lived only in a session scratchpad -- `docs/plan.md:1991`
# recorded the debt ("Productionise it if this is ever re-run"). This is that file,
# repointed and gated. The eta-grid version is NOT reproducible from `scripts/`.
#
# ⚠⚠ THIS RUN IS ONE-BUILD; THE 2026-10-26 RUN WAS NOT. That run read its template VCF,
# cohort array and subjects from `external/PanMixer/starting_data/chr21`, i.e. the
# 340,824-site vg-1.36 tree that `docs/data.md:150` marks STALE, and its releases from a
# 340,824 mirror. Everything here reads OUR 340,849 build, like every other 0.50-0.80
# result. ⚠ CONSEQUENCE: the two read-mapping tables are NOT poolable and the old
# `eta_tau` numbers must not be quoted beside these.
#
# ⚠ Donors are OUR five (HG00138/HG00635/HG01112/HG02698/NA18853), from the paper's
# Methods. Their constants.py runs HG01600 in place of NA18853, and NA18853 (AFR) is
# the LOWEST donor on every metric, so matching their code would move our absolute
# numbers UP and further from their published figures, not closer. Within-pipeline
# deltas are unaffected by the choice.
#
# ⚠ `acc_0p50` IS `eta` = 0 -- no protection, zero target information. It is a SECOND
# control, not a privacy point, so this grid carries TWO baselines and spans SIX
# privacy points, not seven. `docs/data.md:358`.
#
# ⚠ n = 1 RELEASE PER POINT remains the binding limit, exactly as on the eta grid: the
# one movement above noise there (-0.063 pt, all five donors together at eta 500) could
# not be separated from draw variation. A second release per column would have to be a
# second TARGET, which is a different quantity. Do not quote a single column's movement.
#
# ⚠ `set -o pipefail` is LOAD-BEARING and was MISSING from the eta version: the
# six-stage bcftools chain reports only the LAST command's status, so a failure in
# `annotate`, `view -c 1` or `norm` yields a short-but-valid VCF and the run continues.
#
# ⚠ EVERY `rm` USES `${VAR:?}`. The eta version wrote `rm -f $OUT/$G.final.vcf.gz*`
# unquoted and unguarded: with $OUT or $G empty that is `rm -f /*`. `:?` makes the shell
# abort instead of expanding to a top-level path.
set -u -o pipefail
E=/home/ds85/miniconda3/envs/panmixer/bin
export PATH=$E:$PATH
PY=/home/ds85/miniconda3/envs/ranpanmixer/bin/python
P=/data/ds85/RanPanMixer/external/PanMixer
VG=$P/downloaded_tools/vg
REFDIR=$P/starting_data/references
# OUR build -- the single decision that makes this run one-build.
SD=/home/ds85/pm_ourbuild/starting_data/chr21
MIR=/home/ds85/beagle_acc_ourbuild/root
# ⚠ NOT /data: it is 98% full (90 G of 3.7T). Transient footprint is ~1.1 G per graph.
OUT=${OUT:-/home/ds85/readmap_acc}
REPO=/home/ds85/projects/Pangenome/RanPanMixer
D=$REPO/logs/B_readmap_acc
mkdir -p "${OUT:?}" "${D:?}"
DONORS="HG00138 HG00635 HG01112 HG02698 NA18853"
TARGET=HG00438
T_SITES=340849          # OUR axis. Asserted, not assumed -- see GUARD 1.
step(){ echo; echo "######## $* :: $(date +%H:%M:%S)"; }

# graph label -> release row id in the our-build mirror ("base" = unobfuscated)
GRAPHS="base:- acc_0p50:0 acc_0p55:1 acc_0p60:2 acc_0p65:3 acc_0p70:4 acc_0p75:5 acc_0p80:6"
NGRAPH=$(echo $GRAPHS | wc -w)
NDONOR=$(echo $DONORS | wc -w)
NCELL=$((NGRAPH * NDONOR))

step "GUARD 0. the FASTQs, BEFORE the graph loop"
# ⚠⚠ The eta version had `[ -s "$FQ" ] || cat ..._1 ..._2 > "$FQ"` INSIDE the donor loop,
#    which on a reaped FASTQ would WRITE INTO $P/read_fastqs -- the read-only pin. Assert
#    instead, so a missing read set is a clean abort and never an illegal write.
for S in $DONORS; do
  FQ=$P/read_fastqs/${S}.fastq.gz
  [ -s "$FQ" ] || { echo "ABORT: $FQ missing or empty. The eta version would have"; \
                    echo "       WRITTEN it into the pinned tree; this refuses instead."; exit 1; }
  echo "    $S $(stat -Lc%s "$FQ") bytes"
done

step "GUARD 1+2. the axis and the sample order"
echo "grch38#chr21 chr21" > "${OUT:?}/rename_chrs.txt"
$PY - "$T_SITES" "$SD" "$MIR" "$TARGET" <<'PYEOF' || exit 1
import numpy as np, sys, json, pathlib
T, sd, mir, target = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
coh = np.load(f"{sd}/pangenome.npy", mmap_mode="r")
pos = np.load(f"{sd}/pangenome_positions.npy", mmap_mode="r")
subj = np.load(f"{sd}/pangenome_subjects.npy", allow_pickle=True)
# ⚠⚠ GUARD 1. The shape check the eta version had (release vs cohort row) catches a
#    CROSS-axis mix, but NOT the error that actually happened: a run entirely on the
#    stale 340,824 tree, which is internally consistent and silently wrong. Only an
#    assertion that NAMES our axis catches that.
assert coh.shape == (44, T, 2), f"cohort {coh.shape} is not (44, {T}, 2) -- wrong build?"
assert pos.shape[0] == T, f"positions {pos.shape[0]} != {T}"
man = json.loads((pathlib.Path(mir) / "manifest.json").read_text())
assert man["n_sites"] == T, f"mirror n_sites {man['n_sites']} != {T}"
# ⚠ GUARD 2. npy_2_vcf.py:23 rewrites the #CHROM sample names POSITIONALLY from the
#   subjects array, so a template whose columns are in another order would relabel every
#   genotype with no error at all.
assert subj[0] == target, f"subjects[0] is {subj[0]}, not {target}"
rows = [r for r in man["rows"] if r["subject"] == target]
want = ["acc_0p50", "acc_0p55", "acc_0p60", "acc_0p65", "acc_0p70", "acc_0p75", "acc_0p80"]
assert [r["tau"] for r in rows[:7]] == want, \
    f"{target}'s first 7 mirror rows are not the ceiling grid in order: {[r['tau'] for r in rows[:7]]}"
assert [r["row_id"] for r in rows[:7]] == [0, 1, 2, 3, 4, 5, 6], \
    f"{target}'s ceiling rows are not row_id 0..6: {[r['row_id'] for r in rows[:7]]}"
print(f"    axis {coh.shape} OK  positions {pos.shape[0]:,} OK  subjects[0]={subj[0]} OK")
print(f"    mirror n_sites {man['n_sites']:,} OK  {target} rows 0..6 = ceiling grid OK")
PYEOF

step "0. decompress the template VCF once (the converter reads plain text)"
if [ ! -s "${OUT:?}/template.vcf" ]; then
  bcftools view "$SD/pangenome.vcf.gz" -Ov -o "${OUT:?}/template.vcf" || exit 1
fi
TREC=$(grep -vc '^#' "${OUT:?}/template.vcf")
echo "    template: $(stat -c%s "${OUT:?}/template.vcf") bytes, $TREC records"
[ "$TREC" -eq "$T_SITES" ] || { echo "ABORT: template has $TREC records, not $T_SITES"; exit 1; }
# ⚠ GUARD 2b. The template's sample columns must BE the subjects array, in order.
diff <(bcftools query -l "${OUT:?}/template.vcf") \
     <($PY -c "import numpy as np;print('\n'.join(np.load('$SD/pangenome_subjects.npy',allow_pickle=True).tolist()))") \
  > /dev/null || { echo "ABORT: template sample order != pangenome_subjects.npy"; exit 1; }
echo "    template sample order == pangenome_subjects.npy OK"

NREC_FIRST=""
for SPEC in $GRAPHS; do
  G=${SPEC%%:*}; ROW=${SPEC##*:}
  T0=$(date +%s)
  step "GRAPH $G (row ${ROW})"

  # ---- 1. assemble the (44, T, 2) cohort array ------------------------------
  # The converter replaces EVERY genotype column from the npy, so the array must be
  # the whole cohort with only the target's row swapped for our release.
  $PY - "$G" "$ROW" "$SD" "$MIR" "${OUT:?}" "$TARGET" <<'PYEOF' || { echo "ABORT: array assembly failed for $G"; exit 1; }
import numpy as np, sys, pathlib
g, row, sd, mir, out, target = sys.argv[1:7]
out = pathlib.Path(out)
subj = np.load(f"{sd}/pangenome_subjects.npy", allow_pickle=True)
coh = np.array(np.load(f"{sd}/pangenome.npy"))   # (44, T, 2), resident: it is edited
if row != "-":
    rel = np.load(f"{mir}/experiments/exp_1/data/chr21/{row}/new_haplotypes.npy")
    i = int(np.flatnonzero(subj == target)[0])
    if rel.shape != coh[i].shape:
        raise SystemExit(f"shape {rel.shape} != cohort row {coh[i].shape}")
    n_before = int((coh[i] != rel).sum())
    # ⚠⚠ GUARD 3. NOT A NULL EDIT. A flat read-mapping result only means "alignment does
    #    not care" if the graph really was perturbed; otherwise it means "nothing
    #    changed". The eta grid measured 137,359-139,327 of 681,648 (~20%) at every
    #    point, so a floor of 50,000 is loose by 2.7x and still catches an all-missing
    #    release, the wrong target's row, or a mirror row id pointing elsewhere.
    if n_before < 50_000:
        raise SystemExit(f"NULL EDIT: only {n_before:,} of {coh[i].size:,} cells differ "
                         f"from the original. A flat result here would be meaningless.")
    coh[i] = rel
    print(f"    substituted {target} at index {i}; {n_before:,} of "
          f"{coh[i].size:,} cells differ from the original")
else:
    print("    UNOBFUSCATED baseline: cohort array used as-is")
np.save(out / f"{g}.npy", coh)
np.save(out / f"{g}_subjects.npy", subj)
PYEOF

  # ---- 2. npy -> VCF, using THEIR converter unmodified ---------------------
  rm -f "${OUT:?}/template_output.vcf"
  $PY "$P/tools/common/npy_2_vcf.py" "${OUT:?}/${G:?}.npy" "${OUT:?}/template.vcf" \
      > "${D:?}/convert_$G.log" 2>&1
  if [ ! -s "${OUT:?}/template_output.vcf" ]; then
    echo "ABORT: converter produced nothing for $G"; tail -5 "${D:?}/convert_$G.log"; exit 1
  fi
  mv "${OUT:?}/template_output.vcf" "${OUT:?}/${G:?}.raw.vcf"
  rm -f "${OUT:?}/${G:?}.npy" "${OUT:?}/${G:?}_subjects.npy"
  RAW=$(grep -vc '^#' "${OUT:?}/${G:?}.raw.vcf")
  echo "    raw VCF: $RAW records"
  [ "$RAW" -eq "$T_SITES" ] || { echo "ABORT: raw VCF has $RAW records, not $T_SITES"; exit 1; }

  # ---- 3. the bcftools chain, IDENTICAL to scripts/run_vg_baseline.sh ------
  # ⚠ Rename FIRST: `norm -f` needs a contig name a reference FASTA actually provides.
  bcftools annotate --rename-chrs "${OUT:?}/rename_chrs.txt" "${OUT:?}/${G:?}.raw.vcf" -Ou \
   | bcftools view -c 1 -Ou \
   | bcftools norm -N -f "$REFDIR/chr21_ucsc.fa" -Ou \
   | bcftools view --trim-alt-alleles -Ou \
   | bcftools view -t chr21 -v snps -Ou \
   | bcftools sort -T "${OUT:?}/sorttmp" -Oz -o "${OUT:?}/${G:?}.final.vcf.gz" 2>> "${D:?}/prep_$G.log"
  tabix -p vcf -f "${OUT:?}/${G:?}.final.vcf.gz" 2>> "${D:?}/prep_$G.log"
  NREC=$(bcftools index -n "${OUT:?}/${G:?}.final.vcf.gz")
  echo "    graph VCF: $NREC records"
  rm -f "${OUT:?}/${G:?}.raw.vcf"
  # ⚠ GUARD 4. The eta version's floor was `-lt 100000`, loose enough to pass a
  #   half-written VCF. Anchor on the FIRST graph's own count and require every later
  #   graph within 2,000 of it: the eta grid's seven sat in 280,294-280,989, a span of
  #   695. A norm failure or a truncated template moves this by far more.
  if [ -z "$NREC_FIRST" ]; then
    NREC_FIRST=$NREC
    [ "$NREC" -gt 100000 ] || { echo "ABORT: baseline has only $NREC records"; exit 1; }
    echo "    (anchor: later graphs must be within 2,000 of $NREC_FIRST)"
  else
    DIFF=$(( NREC > NREC_FIRST ? NREC - NREC_FIRST : NREC_FIRST - NREC ))
    [ "$DIFF" -le 2000 ] || { echo "ABORT: $NREC records is $DIFF from the anchor $NREC_FIRST"; exit 1; }
  fi

  # ---- 4. index ------------------------------------------------------------
  rm -f "${OUT:?}"/idx*
  $VG autoindex --workflow giraffe --prefix "${OUT:?}/idx" --ref "$REFDIR/chr21_ucsc.fa" \
      --vcf "${OUT:?}/${G:?}.final.vcf.gz" -t 32 > "${D:?}/index_$G.log" 2>&1
  GBZ=${OUT:?}/idx.giraffe.gbz
  [ -s "$GBZ" ] || { echo "ABORT: no index for $G"; tail -5 "${D:?}/index_$G.log"; exit 1; }
  DIST=$(ls "${OUT:?}"/idx*.dist 2>/dev/null | head -1)
  MIN=$(ls "${OUT:?}"/idx*.min 2>/dev/null | head -1)
  ZIP=$(ls "${OUT:?}"/idx*.zipcodes 2>/dev/null | head -1)

  # ---- 5. align the five donors -------------------------------------------
  for S in $DONORS; do
    FQ=$P/read_fastqs/${S}.fastq.gz
    ARGS=(-Z "$GBZ" -d "$DIST" -m "$MIN" -f "$FQ" -t 32)
    [ -n "${ZIP:-}" ] && ARGS+=(-z "$ZIP")
    $VG giraffe "${ARGS[@]}" > "${OUT:?}/tmp.gam" 2>> "${D:?}/align_$G.log"
    $VG stats -a "${OUT:?}/tmp.gam" > "${D:?}/stats_${G}_${S}.txt" 2>> "${D:?}/align_$G.log"
    rm -f "${OUT:?}/tmp.gam"
    A=$(awk -F': ' '/^Total aligned/{print $2}' "${D:?}/stats_${G}_${S}.txt")
    echo "    $S aligned=$A"
  done

  rm -f "${OUT:?}"/idx*
  rm -f "${OUT:?}/${G:?}.final.vcf.gz" "${OUT:?}/${G:?}.final.vcf.gz.tbi"
  echo "    graph $G done in $(( $(date +%s) - T0 ))s"
done

step "COMPLETENESS, counted from disk"
# ⚠ A fan-out that filters its own results cannot tell "found nothing" from "ran
#   nothing". The constant is DERIVED from the grid so the two cannot drift apart --
#   the eta version hardcoded 35 beside a 7-graph list.
N=$(ls "${D:?}"/stats_*.txt 2>/dev/null | wc -l)
echo "stats files: $N (expect $NCELL = $NGRAPH graphs x $NDONOR donors)"
[ "$N" -eq "$NCELL" ] || { echo "INCOMPLETE -- refusing to report"; exit 1; }
rm -f "${OUT:?}/template.vcf"
step "DONE"
