#!/usr/bin/env python
"""ARM 1 in-graph utility: PanMixer's `af_loss` and `ld_loss`, ours and theirs.

⚠⚠ EXPLORATORY, NOT A PRE-REGISTERED ENDPOINT. Decision H (2026-10-08) declared both
NOT adopted, and this was run on 2026-10-14 AFTER outcomes were inspected, at Dylan's
request. Constraint 6 is not violated by computing them -- it is violated by reporting
them as endpoints -- so every table this emits carries the exploratory stamp and they
stay out of the endpoint tables.

⚠⚠ `af_loss` AND `ld_loss` ARE UTILITY METRICS IN PANMIXER'S OWN FRAMING, NOT PRIVACY
METRICS. `af_loss` sits very close to their optimiser's own objective. Their privacy
metrics are `gap_score`, `diploid_gap_score`, `MIA_privacy` and the Beagle reconstruction.
Reporting these as "the privacy metrics they use" would invert their paper.

WHAT THEY MEASURE. Both substitute the released haplotypes into the 44-member pangenome
in place of the target's row and measure the shift against the original pangenome, so
these are IN-GRAPH quantities, not path-level ones:
  af_loss  mean |delta allele frequency| per position, L1, stratified by MAF band
  ld_loss  mean |delta r^2| over SNP pairs inside a sliding window

⚠ THE LD WINDOW IS 5 POSITION UNITS, NOT 5 KB, despite `KB_WINDOW_SIZE = 5`
(ld_loss.py:8) -- it is compared against a raw position difference at ld_loss.py:107.
The median SNP gap on chr21 is 83 bp, so it reaches 2.92% of SNP positions and ~7,450
pairs; a true 5 kb window would be 10,588,565 pairs, about 1,400x more. We run it AS
SHIPPED, because that is the faithful comparison, and state the window every time.

⚠ `af_loss.py:156` reads `row["capacity"]` and never uses it. Our staging writes the
literal "n/a" into that column rather than a number a later reader could mistake for a
real PanMixer setting.

⚠⚠ TWO MIRROR SHAPES, AND THE MANIFEST DECIDES WHICH.
  LEGACY (`logs/B_eta_arm1`, 2026-10-14/25): exp_6 holds EIGHT rows -- HG00438 alone at
  8 capacities -- and `starting_data` is PANMIXER'S vg-1.36 tree, so the comparison is
  cross-build AND single-subject on their side.
  OUR-BUILD (2026-10-28, `manifest['theirs_rows']` present): exp_6 holds 352 rows -- the
  SAME 44 individuals at 8 capacities -- and `starting_data` is OUR 340,849-record tree,
  so both halves are measured against one `pangenome.npy`. Decision F is satisfied on
  the mechanism, and their column finally carries an n.

    python scripts/arm1_af_ld_results.py --root /home/ds85/pm_eval_ourbuild
"""

from __future__ import annotations

import argparse
import json
import math
import statistics as st
from pathlib import Path

ROOT = Path("/data/ds85/RanPanMixer/runs/beagle_attack/root")
TAUS = ("0p00", "0p05", "0p10", "0p25", "0p50", "0p90", "0p99")
TAU_LBL = ["0", "0.05", "0.1", "0.25", "0.5", "0.9", "0.99"]
CAPS = ["0.01", "0.05", "0.1", "0.2", "0.3", "0.5", "0.7", "0.9"]
KEYS = ["all", "complex", "snp_only", "snp_only_maf_0.01", "snp_only_maf_0.05",
        "snp_only_maf_0.1", "snp_only_maf_0.5"]


# ⚠⚠ `wd_af` IS `af_loss` TIMES TWO, AND IT IS THE NUMBER THEIR FIGURES REPORT.
# external/PanMixer/tools/common/gather_results.py:84-90 derives every one of these as
# `2 * af_loss[<stratum>]`; plot_util/utility_in.py:135 then puts the y-label
# "Wasserstein divergence" on them. There is NO Wasserstein or earth-mover computation
# anywhere in their tree -- a grep for wasserstein|earth_mover|emd returns two hits,
# both axis labels. Their Methods Eq. (8) is U_AF = 2*(1/m)*sum_j |p_j - q_j|, and for a
# bi-allelic site the name is exact rather than loose: W1 between two Bernoullis is
# |p-q| and the L1 distance between the two-category vectors is 2|p-q|.
# WE PRINT IT BECAUSE OUR RAW `af_loss` COLUMN IS OFF THEIR PUBLISHED SCALE BY 2x.
# Their key -> our af_loss stratum, exactly as gather_results.py maps them:
WD_MAP = [
    ("wd_af",            "all",               "AF"),
    ("wd_af_snp_only",   "snp_only",          "AF (SNPs only)"),
    ("wd_af_complex",    "complex",           "AF (complex)"),
    ("wd_af_snps_0_01",  "snp_only_maf_0.01", "AF (MAF < 0.01)"),
    ("wd_af_snps_0_05",  "snp_only_maf_0.05", "AF (MAF in [0.01, 0.05))"),
    ("wd_af_snps_0_1",   "snp_only_maf_0.1",  "AF (MAF in [0.05, 0.1))"),
    ("wd_af_snps_0_5",   "snp_only_maf_0.5",  "AF (MAF in [0.1, 0.5])"),
]
# The three their Fig 4A actually plots, with the PUBLISHED value at their
# `eps_private` = 0.002 operating point (docs/plan.md, validation-against-the-paper
# table). ⚠ Panel four of that figure is `ld_loss`, which needs no rescale.
#
# ⚠⚠ THE PUBLISHED COLUMN IS A SCALE REFERENCE, NEVER A MATCHED COMPARISON. These are
# THEIR OWN run of THEIR mechanism on THEIR vg-1.36 build at one chosen privacy level.
# On the our-build mirror their 8-capacity sweep is measured HERE instead, under the
# same `pangenome.npy` as ours, and THAT is the comparison to read. The published
# figures remain printed because they are what a reader of their paper has in mind.
THEIR_FIG4A = {"wd_af": 0.006, "wd_af_snp_only": 0.006, "wd_af_snps_0_05": 0.005}


def label_parts(t: str) -> tuple[str, float]:
    """A release-directory label -> (display string, numeric value).

    ⚠ FOUR NAMESPACES NOW. `acc_0p80` is the 2026-10-27 `acc_ceiling` grid, `eta_250`
    the genome-wide eta grid, `tau_0p99`/`0p99` the two legacy tau forms. The old code
    knew only the last two, so an `acc_` label raised `ValueError: could not convert
    string to float: 'acc_0.80'` inside the label parse -- before any number was read.
    """
    for pre in ("acc_", "eta_", "tau_"):
        if t.startswith(pre):
            v = float(t[len(pre):].replace("p", "."))
            return (f"{v:g}", v)
    v = float(t.replace("p", "."))
    return (f"{v:g}", v)


def _wd_block(cols, get) -> None:
    """Print the `wd_af` family beside the raw `af_loss` it is derived from.

    `cols` is the column keys in order; `get(stratum, col)` returns that cell's raw
    `af_loss`. Rows are emitted into the table already open above.

    ⚠ The MAF strata are DISJOINT BINS, not cumulative thresholds
    (`af_loss.py:42-46`): maf_0.01 is MAF < 0.01, maf_0.05 is [0.01, 0.05),
    maf_0.1 is [0.05, 0.1), maf_0.5 is [0.1, 0.5]. So their own figure label
    "AF (SNPs MAF $<$ 0.05)" on `wd_af_snps_0_05` is WRONG about their own data --
    that panel is the MIDDLE BIN, not every SNP below 0.05, which would be the union
    of two bins. Our labels state the bin.
    """
    for their, ours, figlbl in WD_MAP:
        vals = " | ".join(f"{2.0 * get(ours, c):.6f}" for c in cols)
        pub = THEIR_FIG4A.get(their)
        tail = f" [their Fig 4A at `eps_private`: {pub}]" if pub is not None else ""
        print(f"| `{their}` = 2x [{ours}] -- {figlbl}{tail} | {vals} |")


def _man() -> dict:
    return json.loads((ROOT / "manifest.json").read_text())


def _collect(exp: str, rows: list[dict], point_key: str):
    """One experiment -> (af[stratum][point] = [...], ld[point] = [...], per-subject maps).

    Completeness is the CALLER's business; this reports what is on disk and nothing
    filters silently -- a missing file leaves a short list, which the caller raises on.
    """
    base = ROOT / f"experiments/exp_{exp}/data/chr21"
    pts = []
    for r in rows:
        p = str(r[point_key])
        if p not in pts:
            pts.append(p)
    af = {k: {p: [] for p in pts} for k in KEYS}
    ld = {p: [] for p in pts}
    af_by_sub: dict[str, dict[str, float]] = {}
    ld_by_sub: dict[str, dict[str, float]] = {}
    for r in rows:
        i, p, s = r["row_id"], str(r[point_key]), r["subject"]
        f = base / str(i) / "af_loss.json"
        g = base / str(i) / "ld_loss.txt"
        if f.exists():
            j = json.loads(f.read_text())
            for k in KEYS:
                af[k][p].append(j[k])
            af_by_sub.setdefault(s, {})[p] = j["all"]
        if g.exists():
            a = g.read_text().split()
            v = float(a[0]) / int(a[1])
            ld[p].append(v)
            ld_by_sub.setdefault(s, {})[p] = v
    return af, ld, af_by_sub, ld_by_sub, pts


def _sign_test(diffs: list[float]) -> tuple[int, int, float, float]:
    """(n_negative, n_total, two-sided sign-test p, paired z). Ties dropped, as in
    `scripts/arm1_headtohead.py`, so the two analyses answer the same question."""
    nz = [d for d in diffs if d != 0.0]
    n, k = len(nz), sum(1 for d in nz if d < 0)
    if n == 0:
        return 0, 0, 1.0, 0.0
    # exact two-sided binomial at p = 0.5
    tail = sum(math.comb(n, j) for j in range(min(k, n - k) + 1)) / (2.0 ** n)
    p = min(1.0, 2.0 * tail)
    m, sd = st.mean(nz), (st.stdev(nz) if n > 1 else 0.0)
    z = (m / (sd / math.sqrt(n))) if sd > 0 else 0.0
    return k, n, p, z


def main(argv=None) -> int:
    global ROOT, TAUS, TAU_LBL
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    a = ap.parse_args(argv)
    ROOT = a.root
    man = _man()
    TAUS = tuple(man["taus"])
    TAU_LBL = [label_parts(t)[0] for t in TAUS]

    # OURS is always exp_5, keyed by the release-directory label.
    af, ld, af_by_sub, ld_by_sub, _ = _collect("5", man["rows"], "tau")
    n = [len(af["all"][t]) for t in TAUS]
    if set(n) != {44}:
        raise SystemExit(f"incomplete: n per point = {n}, expected 44 everywhere")

    # THEIRS: 352 rows on the our-build mirror, 8 legacy rows otherwise.
    ourbuild = "theirs_rows" in man
    if ourbuild:
        taf, tld, taf_by_sub, tld_by_sub, tpts = _collect("6", man["theirs_rows"], "capacity")
        tn = [len(taf["all"][p]) for p in tpts]
        if set(tn) != {44}:
            raise SystemExit(f"theirs incomplete: n per capacity = {tn}, expected 44")
        tcaps = tpts
    else:
        base = ROOT / "experiments/exp_6/data/chr21"
        legacy = [json.loads((base / str(i) / "af_loss.json").read_text()) for i in range(8)]
        legacy_ld = [(base / str(i) / "ld_loss.txt").read_text().split() for i in range(8)]
        tcaps = CAPS
        taf = {k: {c: [legacy[i][k]] for i, c in enumerate(CAPS)} for k in KEYS}
        tld = {c: [float(legacy_ld[i][0]) / int(legacy_ld[i][1])] for i, c in enumerate(CAPS)}
        taf_by_sub = tld_by_sub = {}

    build = man.get("build", "PanMixer's vg-1.36 starting_data (CROSS-BUILD)")
    print("ARM 1 IN-GRAPH UTILITY -- EXPLORATORY, NOT PRE-REGISTERED (decision H).")
    print("⚠ `af_loss`/`ld_loss` are UTILITY metrics in PanMixer's framing, NOT privacy.")
    print(f"chr21. build: {build}")
    print(f"Ours: 44 targets x {len(TAUS)} acc_ceiling, median over targets.")
    print(f"Theirs: {'44 subjects' if ourbuild else 'HG00438 only (n=1)'} x "
          f"{len(tcaps)} capacity, median over subjects.")
    print("⚠ ld_loss window is 5 POSITION UNITS, not 5 kb (`ld_loss.py:107`).")
    print("⚠ acc_ceiling 0.50 is eta = 0 -- NO protection; its level is not a privacy result.\n")

    print("**OURS** (columns: `acc_ceiling`; n = 44 targets per cell)\n")
    print("| variable | " + " | ".join(TAU_LBL) + " |")
    print("|---|" + "---|" * len(TAUS))
    for k in KEYS:
        print(f"| `af_loss` [{k}] | " + " | ".join(f"{st.median(af[k][t]):.6f}" for t in TAUS) + " |")
    print("| `ld_loss` per pair | " + " | ".join(f"{st.median(ld[t]):.6f}" for t in TAUS) + " |")
    _wd_block(TAUS, lambda k, t: st.median(af[k][t]))

    print(f"\n**THEIRS** (columns: `capacity`; n = {44 if ourbuild else 1} per cell)\n")
    print("| variable | " + " | ".join(str(c) for c in tcaps) + " |")
    print("|---|" + "---|" * len(tcaps))
    for k in KEYS:
        print(f"| `af_loss` [{k}] | " + " | ".join(f"{st.median(taf[k][c]):.6f}" for c in tcaps) + " |")
    print("| `ld_loss` per pair | " + " | ".join(f"{st.median(tld[c]):.6f}" for c in tcaps) + " |")
    _wd_block(tcaps, lambda k, c: st.median(taf[k][c]))

    # ⚠ THE DIRECTION CHECK IS THE POINT, not the movement size: the pre-registered
    # endpoints move 44/44 the same way. If this one does not, it does not resolve.
    t0, t1 = TAUS[0], TAUS[-1]
    mv = [v[t1] - v[t0] for v in af_by_sub.values() if t0 in v and t1 in v]
    down, up = sum(1 for x in mv if x < 0), sum(1 for x in mv if x > 0)
    print(f"\naf_loss[all] per-target direction, acc_ceiling {TAU_LBL[0]} -> {TAU_LBL[-1]}: "
          f"{down}/{len(mv)} down, {up}/{len(mv)} up, median {st.median(mv):+.3e}")
    print("  -> NO systematic direction, unlike target_fidelity and e_q_u (44/44)."
          if min(down, up) > len(mv) // 4 else "  -> systematic.")

    # ---- the paired head-to-head, available only on the our-build mirror
    if ourbuild and taf_by_sub:
        sat = tcaps[-1]
        print(f"\n**PAIRED BY INDIVIDUAL, ONE BUILD** -- our `acc_ceiling` {TAU_LBL[-1]} "
              f"against their `capacity` {sat}, same 44 people, same `pangenome.npy`.")
        print("⚠⚠ CROSS-OPERATING-POINT: their `capacity` is not a privacy level in our "
              "sense and `acc_ceiling` does not transfer to their mechanism, so this "
              "pairs the two tools at THEIR OWN endpoints, not at matched privacy "
              "(standing constraint 5). It answers 'at each tool's most-perturbing "
              "setting, which disturbs the graph less', nothing more.")
        print("\n| metric | ours (median) | theirs (median) | we win | sign p | paired z |")
        print("|---|---|---|---|---|---|")
        for name, oa, ta in (("`af_loss` [all]", af_by_sub, taf_by_sub),
                             ("`ld_loss` per pair", ld_by_sub, tld_by_sub)):
            subs = sorted(set(oa) & set(ta))
            d = [oa[s][t1] - ta[s][sat] for s in subs if t1 in oa[s] and sat in ta[s]]
            k, nn, p, z = _sign_test(d)
            print(f"| {name} | {st.median([oa[s][t1] for s in subs]):.6f} | "
                  f"{st.median([ta[s][sat] for s in subs]):.6f} | "
                  f"**{k}/{nn}** | {p:.3g} | {z:+.2f} |")
        print("⚠ 'we win' counts individuals where OUR loss is lower. n is the paired "
              "count after dropping exact ties.")

    # ⚠ ON THEIR PUBLISHED SCALE, which is the only scale a reader of their paper has.
    print("\n**ON THEIR PUBLISHED SCALE (`wd_af`), ours against their "
          "`eps_private` point.** ⚠ Cross-operating-point, and the published figures are "
          "from THEIR vg-1.36 build -- see the THEIR_FIG4A note in this script. "
          "A ratio below 1.00 means WE distort less.")
    print("\n| panel | ours (median, worst column) | theirs published | ratio |")
    print("|---|---|---|---|")
    for their, stratum, figlbl in WD_MAP:
        pub = THEIR_FIG4A.get(their)
        if pub is None:
            continue
        worst = max(2.0 * st.median(af[stratum][t]) for t in TAUS)
        print(f"| {figlbl} | {worst:.6f} | {pub:.3f} | **{worst / pub:.2f}x** |")
    print("⚠ 'worst column' is the largest median across the grid, so this is the "
          "least flattering reading of our own sweep.")

    lo = min(st.median(af["all"][t]) for t in TAUS)
    hi = max(st.median(af["all"][t]) for t in TAUS)
    tlo = min(st.median(taf["all"][c]) for c in tcaps)
    thi = max(st.median(taf["all"][c]) for c in tcaps)
    print(f"\naf_loss[all] grid range {lo:.6f}..{hi:.6f} "
          f"(span {hi - lo:.2e}, {100.0 * (hi - lo) / lo:.2f}% of the floor); "
          f"theirs {tlo:.6f}..{thi:.6f}")
    print(f"ld_loss grid range {min(st.median(ld[t]) for t in TAUS):.6f}.."
          f"{max(st.median(ld[t]) for t in TAUS):.6f}; "
          f"theirs {min(st.median(tld[c]) for c in tcaps):.6f}.."
          f"{max(st.median(tld[c]) for c in tcaps):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
