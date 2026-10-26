#!/usr/bin/env python
"""ARM 1 in-graph utility: PanMixer's `af_loss` and `ld_loss`, ours and theirs.

⚠⚠ EXPLORATORY, NOT A PRE-REGISTERED ENDPOINT. Decision H (2026-10-08) declared both
NOT adopted, and this was run on 2026-10-14 AFTER outcomes were inspected, at Dylan's
request. Constraint 6 is not violated by computing them -- it is violated by reporting
them as endpoints -- so every table this emits carries the exploratory stamp and they
stay out of the endpoint tables.

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

    python scripts/arm1_af_ld_results.py
"""

from __future__ import annotations

import argparse
import json
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
# ⚠⚠ THIS IS A CROSS-BUILD, CROSS-OPERATING-POINT COMPARISON AND BOTH HALVES OF THAT
# MATTER. These numbers come from THEIR OWN run of THEIR mechanism on THEIR vg-1.36
# build at one chosen privacy level -- not from the 8-capacity sweep we have staged in
# `exp_6`. The two are visibly not the same object: their published `ld_loss` at
# `eps_private` is 0.003 while `exp_6`'s staged rows measure 0.0380..0.0401, a 13x
# gap, and their published `wd_af` of 0.006 is ABOVE the 0.005468 that 2x their own
# staged maximum (0.002734) gives. So treat the published column as a SCALE REFERENCE
# -- it is what their paper reports and therefore what a reader will have in mind --
# and never as a matched-configuration comparison. Decision F is unsatisfied here.
THEIR_FIG4A = {"wd_af": 0.006, "wd_af_snp_only": 0.006, "wd_af_snps_0_05": 0.005}


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


def _rows():
    return {r["row_id"]: r for r in json.loads((ROOT / "manifest.json").read_text())["rows"]}


def ours():
    base = ROOT / "experiments/exp_5/data/chr21"
    af = {k: {t: [] for t in TAUS} for k in KEYS}
    ld = {t: [] for t in TAUS}
    per_target: dict[str, dict[str, float]] = {}
    for i, r in _rows().items():
        f = base / str(i) / "af_loss.json"
        g = base / str(i) / "ld_loss.txt"
        if f.exists():
            j = json.loads(f.read_text())
            for k in KEYS:
                af[k][r["tau"]].append(j[k])
            per_target.setdefault(r["subject"], {})[r["tau"]] = j["all"]
        if g.exists():
            a = g.read_text().split()
            ld[r["tau"]].append(float(a[0]) / int(a[1]))
    return af, ld, per_target


def theirs():
    base = ROOT / "experiments/exp_6/data/chr21"
    af = [json.loads((base / str(i) / "af_loss.json").read_text()) for i in range(8)]
    ld = [(base / str(i) / "ld_loss.txt").read_text().split() for i in range(8)]
    return af, [float(a[0]) / int(a[1]) for a in ld]


def main(argv=None) -> int:
    # The grid is read from the MIRROR's own manifest, so a second mirror (the
    # exploratory high-`tau` one) is tabled by this same code.
    global ROOT, TAUS, TAU_LBL
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    ROOT = ap.parse_args(argv).root
    TAUS = tuple(json.loads((ROOT / "manifest.json").read_text())["taus"])
    # ⚠ a label may be "eta_250" (the eta namespace) or "0p99" (legacy tau)
    TAU_LBL = [(f"{float(t[4:].replace('p', '.')):g}" if t.startswith("eta_")
                else f"{float(t.replace('p', '.')):g}") for t in TAUS]
    af, ld, per_target = ours()
    n = [len(af["all"][t]) for t in TAUS]
    if set(n) != {44}:
        raise SystemExit(f"incomplete: n per tau = {n}, expected 44 everywhere")
    taf, tld = theirs()

    print("ARM 1 IN-GRAPH UTILITY -- EXPLORATORY, NOT PRE-REGISTERED (decision H).")
    print(f"chr21. Ours: 44 targets x {len(TAUS)} tau, median. Theirs: HG00438 x 8 capacities.")
    print("⚠ ld_loss window is 5 POSITION UNITS, not 5 kb. ⚠ Not same-build.\n")

    print("**OURS**\n")
    print("| variable | " + " | ".join(TAU_LBL) + " |")
    print("|---|" + "---|" * len(TAUS))
    for k in KEYS:
        print(f"| `af_loss` [{k}] | " + " | ".join(f"{st.median(af[k][t]):.6f}" for t in TAUS) + " |")
    print("| `ld_loss` per pair | " + " | ".join(f"{st.median(ld[t]):.6f}" for t in TAUS) + " |")
    _wd_block(TAUS, lambda k, t: st.median(af[k][t]))

    print("\n**THEIRS**\n")
    print("| variable | " + " | ".join(CAPS) + " |")
    print("|---|" + "---|" * 8)
    for k in KEYS:
        print(f"| `af_loss` [{k}] | " + " | ".join(f"{a[k]:.6f}" for a in taf) + " |")
    print("| `ld_loss` per pair | " + " | ".join(f"{x:.6f}" for x in tld) + " |")
    _wd_block(range(len(taf)), lambda k, i: taf[i][k])

    # ⚠ THE DIRECTION CHECK IS THE POINT, not the movement size: the pre-registered
    # endpoints move 44/44 the same way. If this one does not, it does not resolve.
    t0, t1 = TAUS[0], TAUS[-1]
    mv = [v[t1] - v[t0] for v in per_target.values() if t0 in v and t1 in v]
    down, up = sum(1 for x in mv if x < 0), sum(1 for x in mv if x > 0)
    print(f"\naf_loss[all] per-target direction, tau {TAU_LBL[0]} -> {TAU_LBL[-1]}: "
          f"{down}/{len(mv)} down, "
          f"{up}/{len(mv)} up, median {st.median(mv):+.3e}")
    print("  -> NO systematic direction, unlike target_fidelity and e_q_u (44/44)."
          if min(down, up) > len(mv) // 4 else "  -> systematic.")
    # ⚠ ON THEIR PUBLISHED SCALE, which is the only scale a reader of their paper has.
    # Our raw `af_loss` column is off it by 2x, which is why `wd_af` is printed at all.
    print("\n**ON THEIR PUBLISHED SCALE (`wd_af`), ours against their "
          "`eps_private` point.** ⚠ Cross-build and cross-operating-point — see the "
          "THEIR_FIG4A note in this script. A ratio below 1.00 means WE distort less.")
    print("\n| panel | ours (median, worst column) | theirs published | ratio |")
    print("|---|---|---|---|")
    for their, stratum, figlbl in WD_MAP:
        pub = THEIR_FIG4A.get(their)
        if pub is None:
            continue
        worst = max(2.0 * st.median(af[stratum][t]) for t in TAUS)
        print(f"| {figlbl} | {worst:.6f} | {pub:.3f} | "
              f"**{worst / pub:.2f}x** |")
    print("⚠ 'worst column' is the largest median across the grid, so this is the "
          "least flattering reading of our own sweep.")

    lo = min(st.median(af["all"][t]) for t in TAUS)
    hi = max(st.median(af["all"][t]) for t in TAUS)
    print(f"af_loss[all] grid range {lo:.6f}..{hi:.6f}; theirs "
          f"{min(a['all'] for a in taf):.6f}..{max(a['all'] for a in taf):.6f}")
    print(f"ld_loss grid range {min(st.median(ld[t]) for t in TAUS):.6f}.."
          f"{max(st.median(ld[t]) for t in TAUS):.6f}; theirs {min(tld):.6f}..{max(tld):.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
