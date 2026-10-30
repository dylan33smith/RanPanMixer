#!/usr/bin/env python
"""ARM 1's head-to-head FIGURES: the frontier with its spread, the paired difference at
every ceiling, and `wd_af`/`ld_loss` with each tool on its OWN privacy axis.

⚠⚠ WHY THE REFERENCE LINE IS A MEDIAN PLUS A SPREAD AND NOT TWO LINES. Splitting the 44
into "the ones we beat" and "the ones we do not" CONDITIONS ON THE COMPARISON'S OWN
OUTCOME: the two groups are defined by who won, so one line sits below ours and one above
BY CONSTRUCTION, and the figure would prove nothing. The honest forms are (a) one median
with the full spread of all 44 behind it, and (b) the paired difference per individual,
where the 11 we lose appear as a sign change rather than as an invented population. Both
are drawn here. If a split is ever wanted it must key on something fixed BEFORE the
comparison -- their own baseline fidelity, the target's heterozygosity -- never the result.

⚠ `target_fidelity` IS A DIFFER RATE: lower is more faithful. A NEGATIVE paired difference
(ours - theirs) therefore means WE are more faithful. The axis is labelled accordingly
because the sign is the single easiest thing to misread here.

⚠ OURS IS AN EXACT EXPECTATION, THEIRS IS A DRAW. Our `target_fidelity` comes from the
posterior marginals with zero Monte-Carlo error; theirs is read from its own stats.json.
`logs/B_pm_draw_noise` measured their draw sd at ~9.8% of the between-subject spread, so
the residual noise is on their side alone and must be said beside any paired figure.

⚠ THE DENOMINATORS ARE NOT EQUAL and the shares are the quotable numbers, never the raw
counts: ours is K = 86 leave-one-out over its own chain, theirs edits a copy of the
target's own row. That is also why `af_loss` is policy-dependent -- see
`scripts/arm1_afld_sensitivity.py`.

⚠ `capacity` and `acc_ceiling` ARE NOT COMPARABLE PARAMETERS. They never share an x-axis
here. Figure C puts each tool on its own privacy axis for exactly that reason, and the
only way they could share one is a leave-one-out PMI for our releases, which is a separate
measurement and is NOT done here.

⚠ `ld_loss` is parsed exactly as `scripts/arm1_af_ld_results.py:163-164` does
(`float(line0) / int(line1)`), not reimplemented -- two parsers would drift.
"""

from __future__ import annotations

import argparse
import json
import statistics as st
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402

REPO = Path(__file__).resolve().parents[1]
MIRROR = Path("/home/ds85/pm_eval_ourbuild")
OBF = Path("/home/ds85/pm_ourbuild/obf_out")
ENDPOINT = REPO / "results" / "arm1_exact_endpoint.json"
CEILINGS = ("0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80")
# ⚠ "their saturation" is capacity 0.9, the same pairing `af_ld_sensitivity.md` uses.
SATURATION = 0.9
AF_PRIMARY = "snp_only"     # the paper's "SNPs overlapping with 1000 Genomes"
OUT = REPO / "logs" / "B_arm1_figures"


def ours_fidelity() -> dict:
    """target -> {ceiling label: exact target_fidelity}."""
    src = ENDPOINT if ENDPOINT.exists() else Path("/home/ds85/rpm_paper_wt/results/arm1_exact_endpoint.json")
    d = json.loads(src.read_text())
    out = {}
    for rec in d["results"]:
        cells = rec["cells"]
        missing = [c for c in CEILINGS if c not in cells]
        if missing:
            raise SystemExit(f"{rec['target']} is missing ceilings {missing}; a partial "
                             f"frontier would be plotted as if it were whole")
        out[rec["target"]] = {c: float(cells[c]["target_fidelity"]) for c in CEILINGS}
    return out


def theirs_stats() -> dict:
    """subject -> {capacity: {fidelity, n_draws, pmi_gain_frac, ...}}, MEAN over seeds.

    ⚠⚠ FIVE DRAWS, NOT ONE, AND THE CHOICE MOVES THE HEADLINE. Seed 123 lives in
    `obf_out/`, seeds 124-127 in `obf_draws/` (the same split
    `scripts/pm_draw_noise.py:4-9` reads). Against ONE draw we are more faithful in
    **27 of 44** at `acc_ceiling` 0.80; against the 5-draw mean, **33 of 44** -- the
    2026-10-27 entry records that the earlier split verdict was their draw noise
    flipping marginal signs, and that 5 draws settled it. Their within-cell sd at
    saturation is 9.888e-04 against a between-subject sd of 1.007e-02, so draw noise is
    9.8% of the biological spread and the SEM of a 5-draw mean is 4.422e-04.
    ⚠ `changed_frac` is their own field and equals `changed_alleles / total_alleles`
    exactly (verified); it is read rather than recomputed.
    """
    acc: dict = {}
    for pat in (OBF.glob("*/cap*/chr21/stats.json"),
                OBF.with_name("obf_draws").glob("*/cap*/s*/chr21/stats.json")):
        for f in sorted(pat):
            d = json.loads(f.read_text())
            acc.setdefault(d["subject_name"], {}).setdefault(
                round(float(d["capacity"]), 2), {})[int(d["seed"])] = d
    out: dict = {}
    for sub, caps in acc.items():
        for cap, by_seed in caps.items():
            recs = list(by_seed.values())
            out.setdefault(sub, {})[cap] = {
                "fidelity": st.fmean(float(r["changed_frac"]) for r in recs),
                "fidelity_one_draw": float(by_seed[123]["changed_frac"]) if 123 in by_seed else None,
                "n_draws": len(recs),
                "seeds": sorted(by_seed),
                "pmi_gain_frac": st.fmean(float(r["pmi_gain_frac"]) for r in recs),
                "utility_loss_frac": st.fmean(float(r["utility_loss_frac"]) for r in recs),
            }
    return out


def _afld(exp: str, rows: list, point_key: str) -> dict:
    """(subject, point) -> {af: {stratum: v}, ld: v}, read from the per-row outputs."""
    base = MIRROR / "experiments" / exp / "data" / "chr21"
    out: dict = {}
    for r in rows:
        i, sub, pt = r["row_id"], r["subject"], r[point_key]
        rec = {}
        f = base / str(i) / "af_loss.json"
        if f.exists():
            rec["af"] = json.loads(f.read_text())
        g = base / str(i) / "ld_loss.txt"
        if g.exists():
            a = g.read_text().split()
            # ⚠ the SAME two lines as arm1_af_ld_results.py:163-164.
            rec["ld"] = float(a[0]) / int(a[1])
        if rec:
            out[(sub, pt)] = rec
    return out


def afld_both() -> tuple[dict, dict]:
    man = json.loads((MIRROR / "manifest.json").read_text())
    ours = _afld("exp_5", man["rows"], "point_value")
    theirs = _afld("exp_6", man["theirs_rows"], "point_value")
    if not ours or not theirs:
        raise SystemExit(f"af/ld rows missing: ours {len(ours)}, theirs {len(theirs)}")
    return ours, theirs


# ------------------------------------------------------------------ figures

def fig_frontier(ours: dict, theirs: dict, path: Path) -> dict:
    tids = sorted(set(ours) & set(theirs))
    x = [float(c) for c in CEILINGS]
    per = [[ours[t][c] for t in tids] for c in CEILINGS]
    med = [st.median(col) for col in per]
    lo = [min(col) for col in per]
    hi = [max(col) for col in per]
    their_sat = [theirs[t][SATURATION]["fidelity"] for t in tids]

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for t in tids:                                  # every individual, faint
        ax.plot(x, [ours[t][c] for c in CEILINGS], color="#0f5c63", alpha=0.13, lw=0.9)
    ax.fill_between(x, lo, hi, color="#0f5c63", alpha=0.10, lw=0)
    ax.plot(x, med, color="#0f5c63", lw=2.4, marker="o", ms=5,
            label=f"PathGuard, median of {len(tids)}", zorder=5)

    tm = st.median(their_sat)
    ax.axhline(tm, color="#8a5a1b", lw=2.0, ls="--", zorder=4,
               label=f"PanMixer at capacity {SATURATION:g}, median of 5-draw means")
    ax.axhspan(min(their_sat), max(their_sat), color="#8a5a1b", alpha=0.10, lw=0,
               label=f"PanMixer, full spread of {len(tids)}")
    # ⚠ the spread as a rug, so the DISPERSION is visible and not just its envelope
    for v in their_sat:
        ax.plot([x[-1] + 0.006], [v], marker="_", ms=11, color="#8a5a1b", alpha=0.55,
                clip_on=False)

    ax.set_xlabel("acc_ceiling  (PathGuard's knob)")
    ax.set_ylabel("target_fidelity  —  lower is more faithful")
    ax.set_title("Fidelity frontier, chr21, one build, same 44 individuals", fontsize=11)
    ax.set_xticks(x)
    ax.legend(fontsize=8.5, frameon=False, loc="best")
    ax.grid(alpha=0.18, lw=0.6)
    spread = max(their_sat) - min(their_sat)
    effect = abs(med[0] - med[-1])
    fig.text(0.012, 0.055,
             f"⚠ THE SPREAD HERE IS BETWEEN-SUBJECT AND IT SWAMPS THE EFFECT: "
             f"{spread:.4f} across individuals against {effect:.5f} across the whole "
             f"grid, a factor of {spread/effect:.0f}. Unpaired, the two tools are "
             f"indistinguishable at this scale.",
             fontsize=7.4, color="#8c2f28")
    fig.text(0.012, 0.018,
             "The paired figure is what resolves it: same individual, both tools. Ours "
             "is an exact expectation, theirs the mean of 5 draws. acc_ceiling 0.50 is "
             "eta_tau = 0 — a control, not a privacy point.",
             fontsize=7.2, color="#4a5560")
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    fig.savefig(path.with_suffix(".pdf")); fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)
    return {"x": x, "median": med, "min": lo, "max": hi,
            "their_saturation": {t: theirs[t][SATURATION]["fidelity"] for t in tids},
            "their_median": tm}


def fig_paired(ours: dict, theirs: dict, path: Path) -> dict:
    """One panel per ceiling: the per-individual paired difference, sorted."""
    tids = sorted(set(ours) & set(theirs))
    fig, axes = plt.subplots(1, len(CEILINGS), figsize=(16.4, 4.3), sharey=True)
    rec = {}
    for ax, c in zip(axes, CEILINGS):
        d = sorted(ours[t][c] - theirs[t][SATURATION]["fidelity"] for t in tids)
        n_better = sum(1 for v in d if v < 0)
        rec[c] = {"diffs": d, "n_more_faithful": n_better, "n": len(d),
                  "median": st.median(d)}
        cols = ["#0f5c63" if v < 0 else "#8c2f28" for v in d]
        ax.bar(range(len(d)), d, color=cols, width=0.86, linewidth=0)
        ax.axhline(0, color="#161c22", lw=1.0)
        ax.set_title(f"{c}\n{n_better} of {len(d)} more faithful", fontsize=9.5)
        ax.set_xticks([])
        ax.grid(axis="y", alpha=0.18, lw=0.6)
    axes[0].set_ylabel("ours − theirs  (negative = we are more faithful)")
    fig.suptitle("Paired fidelity difference per individual, at every acc_ceiling, "
                 f"against PanMixer at capacity {SATURATION:g}", fontsize=11.5, y=0.99)
    fig.text(0.004, 0.022,
             "One bar per individual, sorted within each panel — bars are NOT the same "
             "person across panels. Teal = PathGuard more faithful. Ours is exact; "
             "theirs is the mean of 5 draws at capacity 0.9.",
             fontsize=7.4, color="#4a5560")
    fig.tight_layout(rect=(0.006, 0.075, 1, 0.94))
    fig.savefig(path.with_suffix(".pdf")); fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)
    return rec


def fig_afld(our_afld: dict, their_afld: dict, theirs: dict, path: Path) -> dict:
    """`wd_af` and `ld_loss`, EACH TOOL ON ITS OWN PRIVACY AXIS (twin x)."""
    o_pts = sorted({p for _, p in our_afld})
    t_pts = sorted({p for _, p in their_afld})
    out = {"ours": {}, "theirs": {}}

    def med_af(src, pt, strat):
        v = [r["af"][strat] for (s, p), r in src.items() if p == pt and "af" in r]
        return st.median(v) if v else None

    def med_ld(src, pt):
        v = [r["ld"] for (s, p), r in src.items() if p == pt and "ld" in r]
        return st.median(v) if v else None

    # ⚠ wd_af = 2 * af_loss IS their published Wasserstein scale (docs/terms.md).
    o_wd = [2 * med_af(our_afld, p, AF_PRIMARY) for p in o_pts]
    t_wd = [2 * med_af(their_afld, p, AF_PRIMARY) for p in t_pts]
    o_ld = [med_ld(our_afld, p) for p in o_pts]
    t_ld = [med_ld(their_afld, p) for p in t_pts]
    # their privacy axis, from their own stats.json
    t_pmi = [st.median([theirs[s][p]["pmi_gain_frac"] for s in theirs if p in theirs[s]])
             for p in t_pts]

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.6))
    for ax, (oy, ty, lab) in zip(axes, ((o_wd, t_wd, "wd_af  (Wasserstein divergence)"),
                                        (o_ld, t_ld, "ld_loss"))):
        ax.plot(o_pts, oy, color="#0f5c63", lw=2.2, marker="o", ms=5, label="PathGuard")
        ax.set_xlabel("acc_ceiling  (PathGuard)", color="#0f5c63")
        ax.tick_params(axis="x", colors="#0f5c63")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.18, lw=0.6)
        top = ax.twiny()
        top.plot(t_pmi, ty, color="#8a5a1b", lw=2.2, marker="s", ms=4.6, ls="--",
                 label="PanMixer")
        top.set_xlabel("pmi_gain_frac  (PanMixer's own privacy axis)", color="#8a5a1b")
        top.tick_params(axis="x", colors="#8a5a1b")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = top.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8.5, frameon=True, framealpha=0.92,
                  edgecolor="none", loc="center left")

    fig.suptitle("In-graph utility, each tool on its OWN privacy axis  —  "
                 f"stratum {AF_PRIMARY}, median over 44, one build", fontsize=11)
    fig.text(0.005, 0.052,
             "⚠ The two x-axes are NOT a shared scale: acc_ceiling bounds "
             "distinguishability, pmi_gain does not, and neither maps to the other. "
             "Reading across them is the comparison this figure refuses to make.",
             fontsize=7.2, color="#8c2f28")
    fig.text(0.005, 0.016,
             "⚠ ld_loss here uses THEIR SHIPPED denominator, under which the gap at "
             "their saturation is ~8x. Under a consistent denominator it is 1.30x — "
             "smaller, but 44/44 rather than 40/44, and the more robust number. Both "
             "belong in the paper. ⚠ af_loss is likewise policy-dependent. "
             "EXPLORATORY under decision H.",
             fontsize=7.2, color="#8c2f28")
    fig.tight_layout(rect=(0, 0.095, 1, 0.94))
    fig.savefig(path.with_suffix(".pdf")); fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)
    out["ours"] = {"points": o_pts, "wd_af": o_wd, "ld_loss": o_ld}
    out["theirs"] = {"capacity": t_pts, "pmi_gain_frac": t_pmi, "wd_af": t_wd,
                     "ld_loss": t_ld}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)

    o = ours_fidelity()
    t = theirs_stats()
    shared = sorted(set(o) & set(t))
    # ⚠ A fan-out that filters its own inputs cannot tell "found nothing" from "ran
    #   nothing". Assert the pairing is complete before anything is plotted.
    if len(shared) != 44:
        raise SystemExit(f"only {len(shared)} individuals appear in BOTH tools "
                         f"(ours {len(o)}, theirs {len(t)}); a partial pairing would be "
                         f"plotted as if it were the cohort")
    nd = set()
    for s in shared:
        if SATURATION not in t[s]:
            raise SystemExit(f"{s} has no capacity {SATURATION:g}; the reference line "
                             f"would be drawn from a different set than the curve")
        nd.add(t[s][SATURATION]["n_draws"])
    # ⚠ An uneven draw count would average different amounts of noise per individual and
    #   the paired signs would shift for a reason that is not the mechanism.
    if len(nd) != 1:
        raise SystemExit(f"uneven draws per individual at saturation: {sorted(nd)}")
    print(f"their reference is the mean of {nd.pop()} draws per cell")
    print(f"44 individuals paired in both tools; ceilings {list(CEILINGS)}; "
          f"their saturation = capacity {SATURATION:g}")

    oa, ta = afld_both()
    print(f"af/ld rows read: ours {len(oa)}, theirs {len(ta)}")

    rec = {"arm": "arm1", "scope": "chr21", "n": len(shared), "ceilings": list(CEILINGS),
           "their_saturation_capacity": SATURATION, "af_stratum": AF_PRIMARY,
           "frontier": fig_frontier(o, t, a.out / "fig_frontier"),
           "paired_by_ceiling": fig_paired(o, t, a.out / "fig_paired_by_ceiling"),
           "afld_own_axes": fig_afld(oa, ta, t, a.out / "fig_afld_own_axes")}
    (a.out / "figures.json").write_text(json.dumps(rec, indent=1) + "\n")

    print("\n**PAIRED FIDELITY DIFFERENCE, PER CEILING.** Metrics are ROWS, "
          "`acc_ceiling` is the COLUMN. n = 44 individuals, paired.\n")
    hdr = ["metric"] + list(CEILINGS)
    print("| " + " | ".join(hdr) + " |"); print("|" + "---|" * len(hdr))
    pb = rec["paired_by_ceiling"]
    print("| more faithful than PanMixer, of 44 | " +
          " | ".join(str(pb[c]["n_more_faithful"]) for c in CEILINGS) + " |")
    print("| median difference (ours - theirs) | " +
          " | ".join(f"{pb[c]['median']:+.6f}" for c in CEILINGS) + " |")
    print("\n⚠ NEGATIVE means WE are more faithful: `target_fidelity` is a differ rate.")
    nd0 = t[shared[0]][SATURATION]["n_draws"]
    print(f"⚠ Ours is an EXACT expectation; theirs is the MEAN OF {nd0} DRAWS at "
          f"capacity {SATURATION:g} (SEM 4.422e-04), so the residual sampling noise is "
          f"on their side alone. Against ONE draw the 0.80 column reads 27 of 44, not "
          f"33: their draw noise flips marginal signs, which is why 5 draws were run.")
    print("⚠ The 0.50 column is `eta_tau` = 0 -- zero target information, a control "
          "rather than a privacy point.")
    print(f"\nFigures -> {a.out}/fig_*.pdf and .png;  numbers -> {a.out}/figures.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
