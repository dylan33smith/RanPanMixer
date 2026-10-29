#!/bin/bash
# ARM 2's GENOME-WIDE utility sweep on the reported `acc_ceiling` 0.50-0.80 grid,
# at all 10 targets, each at ITS OWN fitted `eta_tau`.
#
# ⚠⚠ WHY THIS WAITS. The sweep cannot run at a shared `eta_tau` list. On 2026-10-28 the
# per-target chr21 fit measured `c` = kl/eta^2 spreading **30.9274%** across these same
# ten targets (2.690912e-06 NA20294 to 3.523140e-06 HG01323), against arm 1's 35.7331% --
# so the prediction that arm 2 would be uniform because it holds nothing out was REFUTED.
# `logs/B_eta_grid_v2/genome_uniform.json` is HG00096-only, so reusing its seven `eta`
# would label nine of ten columns with a ceiling that is not theirs. This queue therefore
# BLOCKS on the genome-wide per-target fit and solves each target's own grid from it.
#
# ⚠ Standing constraint 7: utility is comparable only at matched `acc_ceiling`, never at
# matched `eta_tau`. That is the whole reason for the per-target solve.
#
# ⚠ `acc_0p50` is `eta_tau` = 0 -- zero target information, so its LEVEL is the untilted
# floor and not a privacy point. The figure is the MOVEMENT across the grid, per the
# free-lunch rule.
#
#   bash scripts/arm2_genome_sweep_queue.sh
set -u -o pipefail
REPO=/home/ds85/projects/Pangenome/RanPanMixer
PY=/home/ds85/miniconda3/envs/ranpanmixer/bin/python
export PYTHONDONTWRITEBYTECODE=1
D=$REPO/logs/B_acc_arm2
FIT=$D/arm2_genome_kl.json
FIT_STATUS=$D/arm2_genome_kl.status
GRID=$D/arm2_genome_grid.json
OUTD=$REPO/logs/B_acc_arm2/genome_sweep
mkdir -p "${OUTD:?}"
# ⚠ Deliberately modest: the privacy-evaluator run is memory-BANDWIDTH bound and sharing
#   a box with it makes both slower than running either alone. 2 x 11 = 22 cores.
NPAR=${NPAR:-2}
WORKERS=${WORKERS:-11}
step(){ echo; echo "######## $* :: $(date +%H:%M:%S)"; }

step "0. wait for the genome-wide per-target fit"
# ⚠ Poll the SENTINEL, never the log: CLAUDE.md's rule is that a run is not finished
#   until its status file says so.
while [ ! -f "$FIT_STATUS" ]; do sleep 60; done
S=$(cat "$FIT_STATUS")
[ "$S" = "0" ] || { echo "ABORT: the fit exited $S; a partial fit would solve a wrong grid"; exit 1; }
[ -s "$FIT" ] || { echo "ABORT: the fit reported 0 but wrote no JSON at $FIT"; exit 1; }
echo "    fit OK: $FIT"

step "1. solve each target's own seven eta from its own c"
$PY - "$FIT" "$GRID" <<'PYEOF' || exit 1
import json, sys, math, pathlib
sys.path.insert(0, "/home/ds85/projects/Pangenome/RanPanMixer/scripts")
sys.path.insert(0, "/home/ds85/projects/Pangenome/RanPanMixer/src")
# ⚠ IMPORTED, never reimplemented: one definition of the ceiling arithmetic.
from arm1_chr21_kl import CEILINGS, ceiling_from_kl, eta_for_ceiling
fit = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert fit["arm"] == "arm2" and fit["scope"] == "genome", f"{fit['arm']}/{fit['scope']}"
assert fit["knob"] == "eta_tau", f"knob is {fit['knob']}, not eta_tau"
by_t = {}
for r in fit["rows"]:
    by_t.setdefault(r["target"], {})[float(r["knob_value"])] = float(r["kl_exact"])
assert by_t, "no rows in the fit"
out = {"arm": "arm2", "scope": "genome", "weights": fit["weights"],
       "k_states": fit["k_states"], "T": fit["T"], "source_fit": sys.argv[1],
       "ceilings": list(CEILINGS), "per_target": {}}
cs = {}
for t, kls in sorted(by_t.items()):
    pts = [(e, k) for e, k in sorted(kls.items()) if e > 0]
    assert pts, f"{t} has no non-zero eta; c is undefined at eta = 0"
    c_each = [k / (e * e) for e, k in pts]
    c = sum(c_each) / len(c_each)
    resid = (max(c_each) / min(c_each) - 1.0) if min(c_each) > 0 else float("inf")
    cs[t] = c
    out["per_target"][t] = {
        "c": c, "c_residual_within_target": resid,
        "kl_at": {f"{e:g}": k for e, k in pts},
        # 0.50 is eta = 0 exactly, by construction, not by solving.
        "grid": {f"{a:.2f}": (0.0 if a <= 0.5 else eta_for_ceiling(a, c))
                 for a in CEILINGS}}
lo_t, hi_t = min(cs, key=cs.get), max(cs, key=cs.get)
spread = cs[hi_t] / cs[lo_t] - 1.0
out["c_spread_across_targets"] = spread
out["c_min"] = {"target": lo_t, "c": cs[lo_t]}
out["c_max"] = {"target": hi_t, "c": cs[hi_t]}
pathlib.Path(sys.argv[2]).write_text(json.dumps(out, indent=1) + "\n")
print(f"    n targets: {len(cs)}")
print(f"    **`c_quadratic` GENOME-WIDE ACROSS {len(cs)} TARGETS: {cs[lo_t]:.6e} "
      f"({lo_t}) to {cs[hi_t]:.6e} ({hi_t}), a {100*spread:.4f}% spread.**")
print(f"      chr21 at share 1 measured 30.9274%; arm 1's chr21 LOO, 35.7331%.")
print(f"    worst within-target fit conditioning: "
      f"{max(v['c_residual_within_target'] for v in out['per_target'].values()):.3e}")
for t in sorted(out["per_target"]):
    g = out["per_target"][t]["grid"]
    print(f"      {t}: eta 0.80 -> {g['0.80']:,.1f}")
print(f"    GRID -> {sys.argv[2]}")
PYEOF

step "2. the sweep, $NPAR targets at a time x $WORKERS workers each"
TIDS=$($PY -c "import json;print(' '.join(sorted(json.load(open('$GRID'))['per_target'])))")
echo "    targets: $TIDS"
started=0
for T in $TIDS; do
  ETAS=$($PY -c "
import json
g=json.load(open('$GRID'))['per_target']['$T']['grid']
print(' '.join(f\"{g[k]:.10g}\" for k in sorted(g, key=float)))")
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 20; done
  ( $PY -B "$REPO/scripts/eta_sweep.py" --genome --target "$T" --weights uniform \
        --with-kl --workers "$WORKERS" --etas $ETAS \
        --json "${OUTD:?}/${T:?}.json" > "${OUTD:?}/${T:?}.log" 2>&1
    echo $? > "${OUTD:?}/${T:?}.status" ) &
  started=$((started+1))
  echo "    launched $T ($started) with etas: $ETAS"
done
wait

step "3. COMPLETENESS, counted from disk"
# ⚠ A fan-out that filters its own results cannot tell "found nothing" from "ran
#   nothing". Count the files, check every sentinel, and refuse to table a partial grid.
NT=$(echo $TIDS | wc -w)
OK=0; BAD=0
for T in $TIDS; do
  f=${OUTD:?}/${T:?}.status
  if [ -f "$f" ] && [ "$(cat "$f")" = "0" ]; then OK=$((OK+1));
  else BAD=$((BAD+1)); echo "  FAILED: $T -> $(cat "$f" 2>/dev/null || echo no-status)"; fi
done
NJ=$(ls "${OUTD:?}"/*.json 2>/dev/null | wc -l)
echo "sentinels: $OK ok, $BAD bad of $NT; JSONs on disk: $NJ (expect $NT)"
[ "$BAD" -eq 0 ] && [ "$NJ" -eq "$NT" ] || { echo "INCOMPLETE -- refusing to report"; exit 1; }

step "4. the table"
$PY - "$GRID" "${OUTD:?}" <<'PYEOF' || exit 1
import json, pathlib, statistics as st, sys
grid = json.loads(pathlib.Path(sys.argv[1]).read_text())
outd = pathlib.Path(sys.argv[2])
tids = sorted(grid["per_target"])
labels = [f"{c:.2f}" for c in grid["ceilings"]]
per = {}
for t in tids:
    d = json.loads((outd / f"{t}.json").read_text())
    agg = d.get("aggregate") or d.get("rows")
    assert agg and len(agg) == len(labels), f"{t}: {len(agg) if agg else 0} points, want {len(labels)}"
    per[t] = agg
print(f"ARM 2 GENOME-WIDE UTILITY on the reported `acc_ceiling` grid. K = "
      f"{grid['k_states']}, 22 autosomes, T = {grid['T']:,}, weights "
      f"{grid['weights']}, n = {len(tids)} targets, median over targets. Metrics are "
      f"ROWS, `acc_ceiling` is the COLUMN.")
print(f"⚠ EACH TARGET AT ITS OWN `eta_tau`, solved from its own `c_quadratic` "
      f"(spread {100*grid['c_spread_across_targets']:.4f}% across targets) -- matched "
      f"on `acc_ceiling`, which constraint 7 requires, never on `eta_tau`.")
print(f"⚠⚠ The 0.50 column is `eta_tau` = 0: ZERO target information, so its LEVEL is "
      f"the untilted floor and NOT a privacy result. The figure is the MOVEMENT.")
print(f"⚠ `acc_ceiling` is a PINSKER bound; the true ceiling is closer to 0.5 "
      f"than every column label.\n")
hdr = ["metric"] + labels
print("| " + " | ".join(hdr) + " |"); print("|" + "---|" * len(hdr))
ROWS = [("`e_q_u`", "e_q_u", ".9f"), ("`e_match_nonref`", "e_match_nonref", ",.1f"),
        ("`e_differ`", "e_differ", ",.1f"), ("`acc_ceiling` *", "acc_ceiling", ".9f"),
        ("`kl_exact`", "kl_exact", ".6e"), ("`eta_tau` median", "eta", ",.1f")]
for lab, key, fmt in ROWS:
    cells = []
    for i in range(len(labels)):
        vals = [per[t][i].get(key) for t in tids]
        vals = [v for v in vals if v is not None]
        cells.append("n/a (not emitted by this run)" if not vals
                     else f"{st.median(vals):{fmt}}")
    print(f"| {lab} | " + " | ".join(cells) + " |")
mv = [st.median([per[t][i]["e_q_u"] for t in tids]) for i in range(len(labels))]
print(f"\n`e_q_u` MOVEMENT across the whole grid: {mv[-1] - mv[0]:+.6e} "
      f"({mv[0]:.9f} at 0.50 to {mv[-1]:.9f} at 0.80).")
print(f"⚠ The LEVEL at 0.50 is reached by an untilted mosaic carrying zero target "
      f"information, so it is never ours to quote; only the movement is.")
print(f"\nJSON -> {outd}/<target>.json, grid -> {sys.argv[1]}")
PYEOF
step "DONE"
