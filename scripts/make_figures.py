#!/usr/bin/env python
"""The seven figures. Every number is a measured value with a named log.

⚠ PROVENANCE LIVES HERE, NOT ON THE FIGURE. Each series below carries a comment naming the
log it came from. The figures themselves are deliberately bare: short title, short axis
labels, point labels, threshold lines, and nothing else. Explanation is the presenter's job.

    python scripts/make_figures.py [--out paper/figures]
"""

from __future__ import annotations

import argparse
import csv
import json
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
    "font.size": 14, "axes.titlesize": 17, "axes.labelsize": 14,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#babdbe", "xtick.color": REF, "ytick.color": REF,
    "axes.labelcolor": DARK, "text.color": DARK, "figure.dpi": 200,
    "legend.fontsize": 13,
})

TAU = [0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99]

# logs/B_exact_kl_stable.md -- ARM 2, genome-wide, T = 22,134,562, exact
ACC_CEILING = [0.500000000, 0.500000000, 0.500010631, 0.500022737,
               0.500048556, 0.500130655, 0.500234679]

# label, mean beta_v, movement in E_Q[u].  logs/B_toy_genome_demo.md |
# logs/B_beta_concentration.md | logs/B_exact_utility.md | logs/B_arm1_results.md
DILUTION = [("toy genome", 0.111, 0.1975),
            ("chr21", 3.27e-06, 3.636777e-06),
            ("genome, Arm 2", 4.5177e-08, 8.345e-08),
            ("genome, Arm 1", 4.5177e-08, 9.6048e-08)]

# the most successful cell of ANY attack we ran, per tau. logs/B_arm2_attack.md +
# logs/B_arm2_a3a4.md. ⚠ includes A4's multiplicity-inflated block maxima.
BEST_ATTACK = [0.5340, 0.5820, 0.5700, 0.5520, 0.5670, 0.5760, 0.5710]

# tau = 0 -> tau = 0.99 in units of each attack's own standard error.
# ⚠ A4's tau = 0 is exactly 0.5 by construction, so its z is an artefact -- flagged red.
ZMOVE = [("A1", -0.36, False), ("A3", -0.36, False), ("A4", +3.18, True),
         ("Gap Score", -1.45, False), ("diploid Gap Score", +0.32, False),
         ("Beagle", -1.16, False)]

# logs/B_arm1_symmetric.md (utility retained) + logs/B_pm_eval_results.md (gap score)
THEIRS = [(0.01, 86.3770, 2453.0), (0.05, 82.7906, 1565.6), (0.10, 78.1504, 1141.1),
          (0.20, 72.8585, -6514.4), (0.30, 72.7477, -6295.9), (0.50, 72.7455, -6297.7),
          (0.70, 72.7455, -6297.7), (0.90, 72.7455, -6297.7)]
OURS_HH = [(0.00, 73.1055, -6877.4), (0.05, 73.4857, -6088.2), (0.10, 73.0568, -7517.0),
           (0.25, 73.0202, -7976.5), (0.50, 73.3704, -8080.0), (0.90, 73.1317, -7491.8),
           (0.99, 73.0596, -7790.3)]

# logs/B_beta_concentration.md -- positions, movement, headroom
REALLOC = {
    "contiguous": [(1, 3.052265e-03, 0.0041644), (10, 1.654607e-02, 0.1250910),
                   (100, 5.792101e-03, 0.1123222), (1000, 4.030703e-03, 0.1475074),
                   (10000, 3.461291e-04, 0.1400100), (100000, 2.070556e-05, None)],
    "spread": [(1, 0.0, None), (10, 2.470895e-03, None), (100, 6.101752e-04, None),
               (1000, 6.686488e-05, None), (10000, 1.086377e-05, None),
               (100000, 5.079640e-06, None)],
}
PROD = (332521, 3.636777e-06, 0.1032968)


def save(fig, out: Path, name: str):
    fig.savefig(out / f"{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {name}")


def fig1(out):
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    x = np.arange(len(TAU))
    perm = [t / 2 for t in TAU]
    ach = [a - 0.5 for a in ACC_CEILING]
    mp = [i for i, v in enumerate(perm) if v > 0]
    ma = [i for i, v in enumerate(ach) if v > 0]
    ax.plot([x[i] for i in mp], [perm[i] for i in mp], "o-", color=REF, lw=2, ms=7,
            label="permitted")
    ax.plot([x[i] for i in ma], [ach[i] for i in ma], "o-", color=ACCENT, lw=2.4, ms=8,
            label="achieved")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels([f"{t:g}" for t in TAU])
    ax.set_xlabel("$\\tau$"); ax.set_ylabel("attacker advantage")
    ax.set_title("Attacker advantage: permitted vs achieved", loc="left", pad=12)
    ax.annotate("0.495", (x[6], perm[6]), xytext=(0, 10), textcoords="offset points",
                ha="center", color=REF)
    ax.annotate("0.000235", (x[6], ach[6]), xytext=(0, -22), textcoords="offset points",
                ha="center", color=ACCENT)
    ax.legend(frameon=False, loc="lower right")
    save(fig, out, "fig1_permitted_vs_achieved")


def fig2(out):
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    ax.plot([1e-8, 1.0], [1e-8, 1.0], "--", color=GREY, lw=1.2, zorder=1)
    for (lab, b, m), c, dx, dy in zip(DILUTION, [GREY, GREY, ACCENT, ACCENT],
                                      (16, 16, 58, 58), (-22, -22, 16, -24)):
        ax.scatter([b], [m], s=130, color=c, zorder=3, edgecolor="white", linewidth=1.5)
        ax.annotate(lab, (b, m), xytext=(dx, dy), textcoords="offset points",
                    fontsize=13, color=c, ha="left",
                    arrowprops=dict(arrowstyle="-", color=c, lw=0.9) if dx > 40 else None)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("mean $\\beta$ per position"); ax.set_ylabel("movement in utility")
    ax.set_title("Utility movement vs budget per position", loc="left", pad=12)
    save(fig, out, "fig2_budget_dilution")


def fig3(out):
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    x = np.arange(len(TAU))
    bound = [(1 + t) / 2 for t in TAU]
    ax.fill_between(x, 0.5, bound, color=REF, alpha=0.10, zorder=1)
    ax.plot(x, bound, "o-", color=REF, lw=2, ms=6, zorder=3, label="bound")
    ax.axhline(0.5, color=GREY, ls=":", lw=1.4, zorder=2)
    ax.plot(x, BEST_ATTACK, "o-", color=ACCENT, lw=2.4, ms=8, zorder=4, label="best attack")
    ax.scatter([1, 2], [BEST_ATTACK[1], BEST_ATTACK[2]], s=240, facecolor="none",
               edgecolor=BAD, linewidth=2.0, zorder=6)
    ax.annotate("0.582", (x[1], BEST_ATTACK[1]), xytext=(0, 18), textcoords="offset points",
                ha="center", color=ACCENT)
    ax.annotate("0.995", (x[6], bound[6]), xytext=(0, 10), textcoords="offset points",
                ha="center", color=REF)
    ax.annotate("0.571", (x[6], BEST_ATTACK[6]), xytext=(0, -24),
                textcoords="offset points", ha="center", color=ACCENT)
    ax.set_xticks(x); ax.set_xticklabels([f"{t:g}" for t in TAU])
    ax.set_xlabel("$\\tau$"); ax.set_ylabel("attacker accuracy")
    ax.set_ylim(0.47, 1.03)
    ax.set_title("Attacker accuracy vs the bound", loc="left", pad=12)
    ax.legend(frameon=False, loc="upper left")
    save(fig, out, "fig3_attacks_vs_bound")


def fig4(out):
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    labs = [z[0] for z in ZMOVE][::-1]
    vals = [z[1] for z in ZMOVE][::-1]
    flag = [z[2] for z in ZMOVE][::-1]
    y = np.arange(len(labs))
    ax.axvspan(-2, 2, color=REF, alpha=0.09, zorder=1)
    ax.axvline(0, color=REF, lw=1.6, zorder=2)
    for yi, v, f in zip(y, vals, flag):
        ax.plot([0, v], [yi, yi], color=GREY, lw=2, zorder=3)
        ax.scatter([v], [yi], s=130, color=BAD if f else ACCENT, zorder=4,
                   edgecolor="white", linewidth=1.5)
        ax.annotate(f"{v:+.2f}", (v, yi), xytext=(14 if v > 0 else -14, 0),
                    textcoords="offset points", va="center",
                    ha="left" if v > 0 else "right", color=BAD if f else ACCENT)
    ax.set_yticks(y); ax.set_yticklabels(labs)
    ax.set_xlabel("standard errors"); ax.set_xlim(-4.2, 4.6)
    ax.set_ylim(-0.6, len(labs) - 0.4)
    ax.set_title("Change from $\\tau$ = 0 to $\\tau$ = 0.99", loc="left", pad=12)
    save(fig, out, "fig4_tau_zero_control")


def fig5(out):
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    tx = [t[1] for t in THEIRS]; ty = [t[2] for t in THEIRS]
    ox = [o[1] for o in OURS_HH]; oy = [o[2] for o in OURS_HH]
    ax.axhline(0, color=BAD, lw=1.8, zorder=2)
    ax.annotate("identified", (86.9, 0), xytext=(0, 9), textcoords="offset points",
                ha="right", color=BAD, fontsize=13)
    ax.plot(tx, ty, "o-", color=REF, lw=2, ms=8, zorder=3, label="PanMixer")
    for cap, xx, yy, dx, dy in ((0.01, *THEIRS[0][1:], 0, 14), (0.10, *THEIRS[2][1:], 0, 14),
                                (0.20, *THEIRS[3][1:], 16, -4), (0.90, *THEIRS[7][1:], 16, 10)):
        ax.annotate(f"{cap:g}", (xx, yy), xytext=(dx, dy), textcoords="offset points",
                    ha="center" if dx == 0 else "left", fontsize=12, color=REF)
    ax.scatter(ox, oy, s=120, color=ACCENT, zorder=4, edgecolor="white", linewidth=1.5,
               label="RanPanMixer")
    ax.set_ylim(-9200, 3600)
    ax.set_xlabel("utility retained (%)"); ax.set_ylabel("Gap Score")
    ax.set_title("Utility vs identifiability", loc="left", pad=12)
    ax.legend(frameon=False, loc="center right")
    save(fig, out, "fig5_frontier_vs_point")


def fig6(out):
    d = json.loads(FIGDATA.read_text())
    cur = list(csv.DictReader((PUB / "beagle_S2_digitized.csv").open()))
    cx = [float(c["n_obfuscated_variants"]) for c in cur]
    cy = [float(c["accuracy_difference_beagle_minus_majorallele"]) for c in cur]
    lo = [float(c["band_low"]) for c in cur]; hi = [float(c["band_high"]) for c in cur]
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    ax.axhline(0, color=BAD, lw=1.8, zorder=2)
    ax.annotate("baseline", (2000, 0), xytext=(0, 9), textcoords="offset points",
                color=BAD, fontsize=13)
    ax.fill_between(cx, lo, hi, color=REF, alpha=0.16, zorder=1)
    ax.plot(cx, cy, "-", color=REF, lw=2.2, zorder=3, label="PanMixer")
    ox = d["cells_snp"]; oy = d["beagle_diff"]
    ax.scatter(ox, oy, s=120, color=ACCENT, zorder=4, edgecolor="white", linewidth=1.5,
               label="RanPanMixer")
    mx = float(np.median(ox)); at = float(np.interp(mx, cx, cy))
    ax.plot([mx, mx], [at, float(np.median(oy))], ":", color=DARK, lw=1.6, zorder=3)
    ax.annotate(f"{np.median(oy):.3f}", (mx, float(np.median(oy))), xytext=(12, 2),
                textcoords="offset points", color=ACCENT)
    ax.annotate(f"{at:.3f}", (mx, at), xytext=(12, -8), textcoords="offset points",
                color=REF)
    ax.set_xscale("log")
    ax.set_xlabel("number of altered variants")
    ax.set_ylabel("accuracy $-$ baseline")
    ax.set_title("Beagle reconstruction vs baseline", loc="left", pad=12)
    ax.legend(frameon=False, loc="lower left")
    save(fig, out, "fig6_beagle_vs_baseline")


def fig7(out):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.4, 5.4))
    fig.subplots_adjust(wspace=0.28)
    for nm, series in REALLOC.items():
        xs = [s[0] for s in series if s[1] > 0]
        ys = [s[1] for s in series if s[1] > 0]
        ax.plot(xs, ys, "o-", color=ACCENT if nm == "contiguous" else GREY, lw=2.2, ms=7,
                label=nm)
    ax.scatter([PROD[0]], [PROD[1]], s=160, marker="D", color=REF, zorder=5,
               edgecolor="white", linewidth=1.5, label="all positions")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1.2e-06, 9e-02)
    ax.set_xlabel("number of positions"); ax.set_ylabel("movement in utility")
    ax.set_title("Movement", loc="left", pad=12)
    ax.legend(frameon=False, loc="lower left")

    hz = [(s[0], s[1] / s[2]) for s in REALLOC["contiguous"] if s[2]]
    ax2.plot([h[0] for h in hz], [h[1] for h in hz], "o-", color=ACCENT, lw=2.2, ms=7)
    ax2.scatter([PROD[0]], [PROD[1] / PROD[2]], s=160, marker="D", color=REF, zorder=5,
                edgecolor="white", linewidth=1.5)
    ax2.set_xscale("log"); ax2.set_yscale("log"); ax2.set_ylim(1e-05, 4.0)
    ax2.set_xlabel("number of positions"); ax2.set_ylabel("headroom captured")
    ax2.set_title("Headroom captured", loc="left", pad=12)
    save(fig, out, "fig7_budget_reallocation")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=REPO / "paper" / "figures")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    print(f"writing to {a.out}")
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7):
        f(a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
