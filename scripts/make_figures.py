#!/usr/bin/env python
"""The seven figures, from measured values only.

⚠ EVERY NUMBER IN THIS FILE IS A MEASURED VALUE WITH A NAMED LOG. Nothing is illustrative,
nothing is smoothed and nothing is rounded for appearance. The provenance of each series is in
the dict that carries it, and the source line printed on each figure names it.

⚠ Figures 1-2 and 7 are ARM 2 / chr21-scope quantities; 3-6 mix arms and that is stated on the
figure itself, because an unstamped number is uninterpretable rather than merely imprecise
(CLAUDE.md, Results Reporting Format).

    python scripts/make_figures.py [--out paper/figures]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

REPO = Path(__file__).resolve().parents[1]
PUB = Path("/data/ds85/RanPanMixer/runs/beagle_attack/published")
FIGDATA = Path("/data/ds85/RanPanMixer/runs/beagle_attack/figdata.json")

GREY = "#9aa0a6"
DARK = "#202124"
ACCENT = "#1a73e8"
REF = "#5f6368"
BAD = "#c5221f"
plt.rcParams.update({
    "font.size": 12, "axes.titlesize": 14, "axes.labelsize": 12,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#babdbe", "xtick.color": REF, "ytick.color": REF,
    "axes.labelcolor": DARK, "text.color": DARK, "figure.dpi": 200,
})

TAU = [0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99]

# ---- measured series -------------------------------------------------------
# logs/B_exact_kl_stable.md -- ARM 2, genome-wide, 22 autosomes, T = 22,134,562, exact
ACC_CEILING = [0.500000000, 0.500000000, 0.500010631, 0.500022737,
               0.500048556, 0.500130655, 0.500234679]
# logs/B_arm1_results.md -- ARM 1, 44 targets, genome-wide, median, exact
E_DIFFER = [7630132.359, 7630132.278, 7630132.196, 7630131.944,
            7630131.466, 7630129.965, 7630128.055]

# the dilution ladder. mean beta_v against the measured movement in E_Q[u].
DILUTION = [
    # label, mean beta_v, movement, provenance
    ("toy genome\n(2 chromosomes)", 0.111, 0.1975, "logs/B_toy_genome_demo.md"),
    ("chr21, whole budget\non one chromosome", 3.27e-06, 3.636777e-06,
     "logs/B_beta_concentration.md"),
    ("genome-wide\nArm 2, n = 1", 4.5177e-08, 8.345e-08, "logs/B_exact_utility.md"),
    ("genome-wide\nArm 1, n = 44", 4.5177e-08, 9.6048e-08, "logs/B_arm1_results.md"),
]

# the single best (most successful) accuracy achieved by ANY attacker we ran, per tau.
# ⚠ this is the ANTI-US choice: it includes A4's multiplicity-inflated block maxima.
BEST_ATTACK = [0.5340, 0.5820, 0.5700, 0.5520, 0.5670, 0.5760, 0.5710]

# tau = 0 -> tau = 0.99, in units of each attack's OWN standard error.
# ⚠ A4's tau = 0 is EXACTLY 0.5 by construction (eta = 0 makes every block statistic
# identically zero), so its z compares a multiplicity-inflated maximum against an exact
# value. It is shown, not hidden, and annotated.
ZMOVE = [
    ("A1  optimal distinguisher", -0.36, False),
    ("A3  rarity-weighted linkage", -0.36, False),
    ("A4  best block of 650", +3.18, True),
    ("PanMixer  Gap Score", -1.45, False),
    ("PanMixer  diploid Gap Score", +0.32, False),
    ("PanMixer  Beagle vs baseline", -1.16, False),
]

# logs/B_arm1_symmetric.md (utility retained) + logs/B_pm_eval_results.md (gap score).
# chr21, target HG00438, one shared position axis, scored by ONE function.
THEIRS = [  # capacity, % of ceiling retained, genotype_g_to_gstar
    (0.01, 86.3770, 2453.0), (0.05, 82.7906, 1565.6), (0.10, 78.1504, 1141.1),
    (0.20, 72.8585, -6514.4), (0.30, 72.7477, -6295.9), (0.50, 72.7455, -6297.7),
    (0.70, 72.7455, -6297.7), (0.90, 72.7455, -6297.7),
]
OURS_HH = [  # tau, % retained, genotype_g_to_gstar
    (0.00, 73.1055, -6877.4), (0.05, 73.4857, -6088.2), (0.10, 73.0568, -7517.0),
    (0.25, 73.0202, -7976.5), (0.50, 73.3704, -8080.0), (0.90, 73.1317, -7491.8),
    (0.99, 73.0596, -7790.3),
]

# logs/B_beta_concentration.md -- 14 allocations, chr21, beta_share = 1, exact.
# positions carrying the budget -> movement in E_Q[u]; headroom = ceiling - E_Q[u] at tau = 0.
REALLOC = {
    "contiguous": [(1, 3.052265e-03, 0.0041644), (10, 1.654607e-02, 0.1250910),
                   (100, 5.792101e-03, 0.1123222), (1000, 4.030703e-03, 0.1475074),
                   (10000, 3.461291e-04, 0.1400100), (100000, 2.070556e-05, None)],
    "spread out": [(1, 0.0, None), (10, 2.470895e-03, None), (100, 6.101752e-04, None),
                   (1000, 6.686488e-05, None), (10000, 1.086377e-05, None),
                   (100000, 5.079640e-06, None)],
}
PROD = (332521, 3.636777e-06, 0.1032968)   # the production allocation, 1/support on every position


def src(ax, text):
    """The provenance line. ⚠ Drawn with the AXES transform, not data coordinates: the
    first version anchored it at data (0,0), which on a log axis is -inf and silently
    dropped the line from every figure."""
    ax.text(0.0, -0.17, text, transform=ax.transAxes, fontsize=8.5, color=GREY,
            va="top", ha="left", linespacing=1.45)


def save(fig, out: Path, name: str):
    fig.savefig(out / f"{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {name}.png / .pdf")


# ---------------------------------------------------------------- 1
def fig1(out):
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    perm = [t / 2 for t in TAU]
    ach = [a - 0.5 for a in ACC_CEILING]
    x = np.arange(len(TAU))
    mp = [i for i, v in enumerate(perm) if v > 0]
    ax.plot([x[i] for i in mp], [perm[i] for i in mp], "o-", color=REF, lw=2, ms=7,
            label="permitted by the proof")
    m = [i for i, v in enumerate(ach) if v > 0]
    ax.plot([x[i] for i in m], [ach[i] for i in m], "o-", color=ACCENT, lw=2.4, ms=8,
            label="achieved by the mechanism (exact)")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels([f"{t:g}" for t in TAU])
    ax.set_xlabel("privacy setting  $\\tau$")
    ax.set_ylabel("attacker advantage over a coin flip\n(log scale)")
    ax.set_title("We are permitted 0.495 of attacker advantage and we use 0.000235\n"
                 "— 2,109$\\times$ less than the proof allows", loc="left", pad=14)
    for i in (5, 6):
        ax.annotate(f"{perm[i]:.3f}", (x[i], perm[i]), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=10, color=REF)
        ax.annotate(f"{ach[i]:.2e}", (x[i], ach[i]), textcoords="offset points",
                    xytext=(0, -16), ha="center", fontsize=10, color=ACCENT)
    ax.annotate("", xy=(x[6], perm[6]), xytext=(x[6], ach[6]),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=1.2))
    ax.annotate("2,109$\\times$", (x[6], math.sqrt(perm[6] * ach[6])),
                xytext=(-46, 0), textcoords="offset points", fontsize=11, color=DARK)
    ax.annotate("$\\tau$ = 0 and 0.05: leakage is below what the\ncomputation resolves, "
                "reported as 0", (x[1], 3e-5), fontsize=9, color=GREY)
    ax.legend(frameon=False, loc="lower right", fontsize=11)
    src(ax, "ARM 2, genome-wide, 22 autosomes, T = 22,134,562, n = 1 target, exact "
            "(no sampling error).  Provenance: logs/B_exact_kl_stable.md\n"
            "Advantage = accuracy − 0.5.  The permitted line is (1+τ)/2 − 0.5 = τ/2, "
            "i.e. Theorem 1's p_succ_bound.")
    save(fig, out, "fig1_permitted_vs_achieved")


# ---------------------------------------------------------------- 2
def fig2(out):
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    b = [d[1] for d in DILUTION]; mv = [d[2] for d in DILUTION]
    lo, hi = 1e-8, 1.0
    ax.plot([lo, hi], [lo, hi], "--", color=GREY, lw=1.2, zorder=1)
    ax.annotate("movement = per-position budget\n(slope 1)", (3e-4, 3e-4), rotation=30,
                fontsize=9, color=GREY, ha="center")
    cols = [GREY, GREY, ACCENT, ACCENT]
    for (lab, bb, mm, _), c in zip(DILUTION, cols):
        ax.scatter([bb], [mm], s=130, color=c, zorder=3, edgecolor="white", linewidth=1.5)
    ax.annotate(DILUTION[0][0], (b[0], mv[0]), xytext=(-12, -42), textcoords="offset points",
                fontsize=10, color=REF, ha="center")
    # ⚠ the two genome-wide points sit almost on top of each other (8.3e-08 against
    # 9.6e-08), so their labels are STACKED below with leaders rather than offset from
    # each point independently, which drew them over the chr21 label.
    ax.annotate(DILUTION[1][0], (b[1], mv[1]), xytext=(12, 16), textcoords="offset points",
                fontsize=10, color=REF)
    ax.annotate("genome-wide, Arm 2 (n = 1): 8.345e−08", (b[2], mv[2]), xytext=(58, 18),
                textcoords="offset points", fontsize=10, color=ACCENT,
                arrowprops=dict(arrowstyle="-", color=ACCENT, lw=0.9))
    ax.annotate("genome-wide, Arm 1 (n = 44): 9.605e−08", (b[3], mv[3]), xytext=(58, -26),
                textcoords="offset points", fontsize=10, color=ACCENT,
                arrowprops=dict(arrowstyle="-", color=ACCENT, lw=0.9))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("mean privacy budget per genome position,  $\\beta_v$")
    ax.set_ylabel("measured movement in expected utility\nacross the whole $\\tau$ range")
    ax.set_title("The effect is budget dilution: movement tracks the per-position budget\n"
                 "across six orders of magnitude — and the two arms land together",
                 loc="left", pad=14)
    src(ax, "Each point is an exact measurement on a different scope. The two genome-wide "
            "points are INDEPENDENT experiments — a single external target (Arm 2) and 44 "
            "held-out cohort members (Arm 1) — measured on\ndifferent days by different code "
            "paths, agreeing to 15%.  Provenance: logs/B_toy_genome_demo.md, "
            "logs/B_beta_concentration.md, logs/B_exact_utility.md, logs/B_arm1_results.md")
    save(fig, out, "fig2_budget_dilution")


# ---------------------------------------------------------------- 3
def fig3(out):
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    x = np.arange(len(TAU))
    bound = [(1 + t) / 2 for t in TAU]
    ax.fill_between(x, 0.5, bound, color=REF, alpha=0.10, zorder=1)
    ax.plot(x, bound, "o-", color=REF, lw=2, ms=6, zorder=3,
            label="p_succ_bound — what the proof permits")
    ax.axhline(0.5, color=GREY, ls=":", lw=1.4, zorder=2)
    ax.annotate("a coin flip, 0.5", (0.05, 0.5), xytext=(0, 7), textcoords="offset points",
                fontsize=9.5, color=GREY)
    ax.plot(x, BEST_ATTACK, "o-", color=ACCENT, lw=2.4, ms=8, zorder=4,
            label="best result of ANY of our 8 attacks")
    for i in (1,):
        ax.annotate(f"{BEST_ATTACK[i]:.4f}", (x[i], BEST_ATTACK[i]),
                    xytext=(0, 11), textcoords="offset points", ha="center",
                    fontsize=10.5, color=ACCENT)
    ax.annotate(f"{bound[6]:.4f}", (x[6], bound[6]), xytext=(0, 9),
                textcoords="offset points", ha="center", fontsize=10.5, color=REF)
    ax.annotate(f"{BEST_ATTACK[6]:.4f}", (x[6], BEST_ATTACK[6]), xytext=(0, -20),
                textcoords="offset points", ha="center", fontsize=10.5, color=ACCENT)
    ax.annotate("permitted but unused", (3.2, 0.72), fontsize=11, color=REF, alpha=0.8)
    ax.annotate("the attack is ABOVE the bound here.\nThe experiment does not resolve below it\n"
                "at these two $\\tau$ (minimum detectable\nadvantage 0.0556) \u2014 reported, "
                "not claimed.\nThe PROOF carries these two columns.",
                (1.0, 0.5820), xytext=(-22, 150), textcoords="offset points", fontsize=9.5,
                color=BAD, ha="left",
                arrowprops=dict(arrowstyle="-", color=BAD, lw=1))
    ax.scatter([1, 2], [BEST_ATTACK[1], BEST_ATTACK[2]], s=210, facecolor="none",
               edgecolor=BAD, linewidth=1.8, zorder=6)
    ax.set_xticks(x); ax.set_xticklabels([f"{t:g}" for t in TAU])
    ax.set_xlabel("privacy setting  $\\tau$"); ax.set_ylabel("attacker accuracy")
    ax.set_ylim(0.47, 1.02)
    # ⚠ the title must not paper over the two cells where the point estimate EXCEEDS the
    # bound. They are the documented not-significant excesses at tau = 0.05 and 0.1, and a
    # title claiming "no attack comes near the bound" is contradicted by the figure itself.
    ax.set_title("From $\\tau$ = 0.25 the bound runs away from every attack; at $\\tau$ = 0.05 "
                 "and 0.1\nthe best of eight EXCEEDS it by point estimate \u2014 not significantly",
                 loc="left", pad=14)
    ax.legend(frameon=False, loc="upper center", fontsize=11)
    src(ax, "The blue line is the MOST SUCCESSFUL cell of every attack we ran at each τ — "
            "A1, A3 and A4 across three decoy difficulties and three block sizes — which is "
            "the choice least flattering to us.\n⚠ It includes A4's multiplicity-inflated "
            "block maxima: picking the best of 650 blocks is optimistically biased, and after "
            "Bonferroni correction the largest excess is z = +2.65 against a critical value of "
            "4.51.\nARM 2, chr21 at full budget, 250 trials per cell per hypothesis.  "
            "Provenance: logs/B_arm2_attack.md, logs/B_arm2_a3a4.md")
    save(fig, out, "fig3_attacks_vs_bound")


# ---------------------------------------------------------------- 4
def fig4(out):
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    labs = [z[0] for z in ZMOVE][::-1]
    vals = [z[1] for z in ZMOVE][::-1]
    flag = [z[2] for z in ZMOVE][::-1]
    y = np.arange(len(labs))
    ax.axvspan(-2, 2, color=REF, alpha=0.09, zorder=1)
    ax.axvline(0, color=REF, lw=1.6, zorder=2)
    for yi, v, f in zip(y, vals, flag):
        ax.plot([0, v], [yi, yi], color=GREY, lw=2, zorder=3)
        ax.scatter([v], [yi], s=120, color=BAD if f else ACCENT, zorder=4,
                   edgecolor="white", linewidth=1.5)
        ax.annotate(f"{v:+.2f}", (v, yi), xytext=(14 if v > 0 else -14, 0),
                    textcoords="offset points", va="center",
                    ha="left" if v > 0 else "right", fontsize=10.5,
                    color=BAD if f else ACCENT)
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=11)
    ax.set_xlabel("change in the attack's outcome from $\\tau$ = 0 to $\\tau$ = 0.99,\n"
                  "in units of that attack's own standard error")
    ax.set_xlim(-4.2, 4.6); ax.set_ylim(-0.7, len(labs) - 0.3)
    ax.annotate("±2 standard errors", (2.0, -0.55), xytext=(6, 0),
                textcoords="offset points", fontsize=9.5, color=REF)
    ax.set_title("Turning the mechanism from OFF to MAXIMUM moves five of six attacks by\n"
                 "less than 1.5 standard errors — these nulls do not test the tilt",
                 loc="left", pad=14)
    src(ax, "τ = 0 applies NO protection: the release is an untilted mosaic of other "
            "people's haplotypes carrying zero target information. An attack that fails "
            "equally at τ = 0 and τ = 0.99 is therefore failing\nagainst the release "
            "FORM, not against the privacy mechanism. Standing constraint 7: a failed attack "
            "is a LOWER bound on distinguishability; only the proof gives an upper bound "
            "(figure 1).\n⚠ RED: A4's τ = 0 value is EXACTLY 0.5 by construction "
            "(η = 0 makes every block statistic identically zero), so its z compares a "
            "multiplicity-inflated maximum against an exact value rather than two\nsampled "
            "ones. Shown rather than dropped. Its proper test is the corrected bound check: "
            "z = +2.65 against a critical value of 4.51.")
    save(fig, out, "fig4_tau_zero_control")


# ---------------------------------------------------------------- 5
def fig5(out):
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    tx = [t[1] for t in THEIRS]; ty = [t[2] for t in THEIRS]
    ox = [o[1] for o in OURS_HH]; oy = [o[2] for o in OURS_HH]
    ax.axhline(0, color=BAD, lw=1.8, zorder=2)
    ax.annotate("identification threshold — above this line the target is re-identified",
                (72.6, 0), xytext=(0, 8), textcoords="offset points", fontsize=10, color=BAD)
    ax.plot(tx, ty, "o-", color=REF, lw=2, ms=8, zorder=3, label="PanMixer, 8 capacities")
    # ⚠ label offsets are per-point: the 0.2 and 0.9 capacities sit almost on top of each
    # other at the saturated end, and one shared offset drew them over one another.
    for cap, xx, yy, dx, dy in ((0.01, *THEIRS[0][1:], 8, 10), (0.10, *THEIRS[2][1:], 8, 10),
                                (0.20, *THEIRS[3][1:], 14, -4), (0.90, *THEIRS[7][1:], 14, 12)):
        ax.annotate(f"capacity {cap:g}", (xx, yy), xytext=(dx, dy),
                    textcoords="offset points", fontsize=9.5, color=REF)
    ax.scatter(ox, oy, s=120, color=ACCENT, zorder=4, edgecolor="white", linewidth=1.5,
               label="RanPanMixer, 7 $\\tau$ values")
    ax.set_ylim(-9600, 4200)
    ax.annotate("all seven $\\tau$ sit here,\nspanning 0.47 points of utility",
                (float(np.mean(ox)), min(oy)), xytext=(26, -34), textcoords="offset points",
                fontsize=10.5, color=ACCENT,
                arrowprops=dict(arrowstyle="-", color=ACCENT, lw=1))
    ax.annotate("", xy=(86.4, 3500), xytext=(72.7, 3500),
                arrowprops=dict(arrowstyle="<->", color=REF, lw=1.2))
    ax.annotate("their frontier spans 13.63 points of utility", (79.5, 3650), ha="center",
                fontsize=10.5, color=REF)
    ax.set_xlabel("utility retained, as a % of the achievable ceiling on a shared axis")
    ax.set_ylabel("Gap Score  (their own linkage attack)\npositive = target identified")
    ax.set_title("PanMixer can trade utility for identifiability; we cannot reach that trade\n"
                 "— a frontier against a point at its conservative end", loc="left", pad=14)
    ax.legend(frameon=False, loc="center right", fontsize=11)
    src(ax, "chr21, target HG00438, both tools' releases scored by ONE function on ONE shared "
            "position axis of 331,416 positions. Ceiling 0.014951984.\n⚠ NOT a same-build "
            "comparison: their releases come from their own vg-1.36 pipeline and ours from our "
            "vg-1.68 rebuild, and the two builds disagree about the target's own genome at "
            "1.34% of cells —\nwhich is LARGER than the gap between the tools. Only each "
            "side's share of its own total is quoted.  ⚠ Their identification at low "
            "capacity is NOT a defect: the knob is meant to be turned up, and locating that\n"
            "threshold is what their published sweep is for.  Provenance: "
            "logs/B_arm1_symmetric.md, logs/B_pm_eval_results.md")
    save(fig, out, "fig5_frontier_vs_point")


# ---------------------------------------------------------------- 6
def fig6(out):
    d = json.loads(FIGDATA.read_text())
    cur = list(csv.DictReader((PUB / "beagle_S2_digitized.csv").open()))
    cx = [float(c["n_obfuscated_variants"]) for c in cur]
    cy = [float(c["accuracy_difference_beagle_minus_majorallele"]) for c in cur]
    lo = [float(c["band_low"]) for c in cur]; hi = [float(c["band_high"]) for c in cur]
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.axhline(0, color=BAD, lw=1.8, zorder=2)
    ax.annotate("parity with the trivial baseline — above this line the attack beats it",
                (2200, 0), xytext=(0, 8), textcoords="offset points", fontsize=10, color=BAD)
    ax.fill_between(cx, lo, hi, color=REF, alpha=0.16, zorder=1)
    ax.plot(cx, cy, "-", color=REF, lw=2.2, zorder=3, label="PanMixer, their published curve")
    ox = d["cells_snp"]; oy = d["beagle_diff"]
    ax.scatter(ox, oy, s=120, color=ACCENT, zorder=4, edgecolor="white", linewidth=1.5,
               label="RanPanMixer, 7 $\\tau$ values (44 targets each)")
    at = float(np.interp(np.median(ox), cx, cy))
    ax.plot([np.median(ox), np.median(ox)], [at, np.median(oy)], ":", color=DARK, lw=1.6,
            zorder=3)
    ax.annotate(f"at a matched {np.median(ox):,.0f} altered variants:\n"
                f"theirs {at:.3f}, ours {np.median(oy):.3f}",
                (np.median(ox), (at + np.median(oy)) / 2), xytext=(16, -6),
                textcoords="offset points", fontsize=10.5, color=DARK)
    ax.set_xscale("log")
    ax.set_xlabel("number of altered variants  (their figure's own x-axis)")
    ax.set_ylabel("Beagle accuracy $-$ trivial-baseline accuracy\nnegative = the attack loses")
    # ⚠ computed, not typed: the first version hardcoded "0.145" in the title while the
    # annotation computed -0.149 from the curve, so the figure disagreed with itself.
    ax.set_title(f"Beagle loses to a trivial baseline against both tools\n"
                 f"— ours by {abs(np.median(oy)):.3f}, theirs by {abs(at):.3f} at a matched "
                 f"number of altered variants", loc="left", pad=14)
    ax.legend(frameon=False, loc="lower left", fontsize=11)
    src(ax, "The baseline is \"guess the commonest 1000 Genomes allele at every site\", which "
            "is PanMixer's own published comparator. Subset snp_only (258,610 sites), the "
            "paper's \"SNPs overlapping with 1000 Genomes\".\n⚠ THEIR CURVE IS OUR "
            "DIGITIZATION OF THEIR FIGURE, NOT A PUBLISHED TABLE: the published paper prints "
            "no Beagle number at all (its Supplementary Fig. S4 is cited and absent), so the "
            "49 points are read from the\nvector path of the PREPRINT's Fig. S2 and were "
            "independently re-extracted with a second tool (y agrees to 7e-05, x to 0.007% of "
            "the axis).  ⚠ Our points are their CODE's protocol, which excludes nobody "
            "from the\nreference panel; their PAPER excludes the target individuals, and that "
            "run is still going.  Provenance: logs/B_arm1_attacks.md, "
            "runs/beagle_attack/published/PROVENANCE.md")
    save(fig, out, "fig6_beagle_vs_baseline")


# ---------------------------------------------------------------- 7
def fig7(out):
    # ⚠ laid out explicitly. The first version let two panel titles and three annotations
    # collide: a figure whose words sit on top of each other is unreadable regardless of
    # whether its numbers are right.
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14.6, 6.2))
    fig.subplots_adjust(wspace=0.30, top=0.80)
    for nm, series in REALLOC.items():
        xs = [s_[0] for s_ in series if s_[1] > 0]
        ys = [s_[1] for s_ in series if s_[1] > 0]
        col = ACCENT if nm == "contiguous" else GREY
        ax.plot(xs, ys, "o-", color=col, lw=2.2, ms=7, label=f"{nm} positions")
    # the production allocation goes in the LEGEND rather than as a floating label: its
    # callout collided with the spread-series annotation in the previous layout.
    ax.scatter([PROD[0]], [PROD[1]], s=160, marker="D", color=REF, zorder=5,
               edgecolor="white", linewidth=1.5,
               label="production allocation\n(all 332,521 positions, 1/support)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(1.2e-06, 9e-02)
    ax.annotate("10 contiguous positions: 1.65e\u221202\n4,550$\\times$ the production allocation",
                (10, 1.654607e-02), xytext=(24, -6), textcoords="offset points",
                fontsize=10.5, color=ACCENT)
    ax.annotate("1 spread-out position moves NOTHING:\nthat index is a no-call for this\n"
                "target \u2014 the honest cost of choosing\npositions without looking at it",
                (10, 2.470895e-03), xytext=(-2, -96), textcoords="offset points",
                fontsize=9.5, color=GREY, ha="left",
                arrowprops=dict(arrowstyle="-", color=GREY, lw=1))
    ax.set_xlabel("positions carrying the whole privacy budget")
    ax.set_ylabel("movement in expected utility\nacross the $\\tau$ range")
    ax.set_title("Concentrating the budget buys up to\n4,550$\\times$ the movement — and is privacy-NEUTRAL",
                 loc="left", fontsize=13, pad=10)
    ax.legend(frameon=False, loc="lower right", fontsize=10, labelspacing=0.9)

    hz = [(s_[0], s_[1] / s_[2]) for s_ in REALLOC["contiguous"] if s_[2]]
    ax2.plot([h[0] for h in hz], [h[1] for h in hz], "o-", color=ACCENT, lw=2.2, ms=7)
    ax2.scatter([PROD[0]], [PROD[1] / PROD[2]], s=160, marker="D", color=REF, zorder=5,
                edgecolor="white", linewidth=1.5)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_ylim(1e-05, 4.0)
    ax2.annotate("1 position captures 73% of its headroom —\nbut that headroom is only 0.0042, because\n"
                 "the baseline already matches that position\n99.6% of the time",
                 (1, 0.733), xytext=(30, -78), textcoords="offset points", fontsize=9.5,
                 color=DARK, arrowprops=dict(arrowstyle="-", color=DARK, lw=1))
    ax2.annotate("production allocation:\n0.0035% of its headroom", (PROD[0], PROD[1] / PROD[2]),
                 xytext=(-14, 34), textcoords="offset points", fontsize=10, color=REF,
                 ha="right", arrowprops=dict(arrowstyle="-", color=REF, lw=1))
    ax2.set_xlabel("positions carrying the whole privacy budget")
    ax2.set_ylabel("share of the AVAILABLE headroom captured")
    ax2.set_title("But the positions where the tilt has power\nare the ones that needed it least",
                  loc="left", fontsize=13, pad=10)
    src(ax, "chr21, T = 332,521, K = 88, target HG00096, ARM 2, whole budget on one "
            "chromosome. 14 allocations \u00d7 3 \u03c4, exact.\n\u26a0 EVERY allocation sums "
            "to 1, so \u03b4_u = 1 and \u03c4, \u03b7_\u03c4, \u03b5_\u03c4 and "
            "p_succ_bound are IDENTICAL on every row: reallocation is privacy-neutral.\n"
            "\u26a0 It does not improve utility, it REDEFINES it \u2014 each allocation is a "
            "different utility function, so levels are not comparable across rows "
            "(constraint 8).\n\u26a0 The positions must be chosen WITHOUT looking at the "
            "target, or Theorem 1 is void (constraint 1).  \u26a0 \u03b2 is a pre-registered "
            "field: changing it needs a NEW dated pre-registration.  "
            "Provenance: logs/B_beta_concentration.md")
    save(fig, out, "fig7_budget_reallocation")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=REPO / "paper" / "figures")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    print(f"writing to {a.out}")
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7):
        f(a.out)
    print(f"\n7 figures, PNG + PDF, in {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
