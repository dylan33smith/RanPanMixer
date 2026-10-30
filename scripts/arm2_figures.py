#!/usr/bin/env python3
"""ARM 2's measured attacker accuracy against the acc line, with error bars.

WHY THIS FILE EXISTS. The paper's only accuracy-vs-bound figure was
`fig3_attacks_vs_bound`, drawn on the dead `tau` grid against `p_succ_bound`. That
bound is 1.0 to float64 at the `eta_tau` values the 0.50-0.80 grid uses, so the old
figure cannot be relabelled onto the current regime -- it has to be rebuilt. This is
that rebuild, and it is the figure Dylan's 2026-10-29 rule is about.

THE LINE IS THE DIAGONAL, not a horizontal. In ceiling mode each target runs at its
own `eta_tau`, solved so that its measured `acc_ceiling` lands on the column label.
So the bound at column c is c itself and the permitted region is y <= x. That is worth
stating because every earlier figure in this project drew a horizontal bound.

THE RULE THIS FIGURE ENFORCES. Every plotted point carries its ERROR BAR, 1.96 x SE
either side, and where a point sits ABOVE the line its bar must be shown reaching
BELOW it. `_gate_bars_reach_below` asserts exactly that and fails the run otherwise --
a figure that can draw an unqualified exceedance is the thing being prevented.

⚠ ARM 2 ONLY, chr21, share 1.0. Arm 2's utility sweep is genome-wide and the two may
never share a panel (standing constraint 5). Arm 1's ceiling is PER PAIR, so its
reference is a band rather than this diagonal; it needs its own figure.
⚠ The three decoy conditions share each release and must not be pooled, so they are
drawn as three series and never averaged.
⚠ Pinsker is an inequality: the TRUE ceiling is closer to 0.5 than the diagonal, and
the looseness runs in the direction that flatters us. Said on the figure, not implied.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from arm2_attack import Z95, classify_excess                      # noqa: E402

A1_JSON = REPO / "logs" / "B_acc_arm2" / "a1_ceil.json"
A34_JSON = REPO / "logs" / "B_acc_arm2" / "a34_ceil.json"
OUT = REPO / "logs" / "B_acc_arm2" / "figures"

CEILINGS = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80)
CONDS = ("D1_extremal", "D2_nearest", "D3_random")
# A3 over its 21 cells, A4 over its 15,603 block-cells. Read from the run, not assumed.
MARK = {"D1_extremal": "o", "D2_nearest": "s", "D3_random": "^"}
COL = {"D1_extremal": "#1b4965", "D2_nearest": "#bc4749", "D3_random": "#6a994e"}


def _k(d: dict, c: float):
    """A JSON round-trip leaves the column label as a string key."""
    for cand in (f"{c:g}", f"{c:.2f}", c):
        if cand in d:
            return d[cand]
    raise KeyError(f"no column {c!r} in {sorted(d)[:9]}")


def load() -> dict:
    """-> {attack: {cond: [(ceiling, bal, se), ...]}}, plus the corrected thresholds.

    ⚠ FAILS LOUDLY on a missing cell. A figure that silently drops a point cannot
    tell "below the line" from "not measured" (the fan-out rule in CLAUDE.md).
    """
    if not A1_JSON.exists() or not A34_JSON.exists():
        raise SystemExit(f"need both {A1_JSON} and {A34_JSON}")
    a1 = json.loads(A1_JSON.read_text())
    a34 = json.loads(A34_JSON.read_text())
    if a1.get("knob") != "acc_ceiling":
        raise SystemExit(f"{A1_JSON} is not a ceiling run: knob={a1.get('knob')!r}")

    series: dict = {}
    series["A1 — optimal distinguisher"] = {
        cond: [(c, (v := _k(a1["accuracy"][cond], c))["balanced"], v["se_balanced"])
               for c in CEILINGS]
        for cond in CONDS}
    series["A3 — rarity-weighted linkage"] = {
        cond: [(c, (v := _k(a34["table"]["a3"][cond], c))["balanced"], v["se"])
               for c in CEILINGS]
        for cond in CONDS}
    # ⚠ A4 is a MAX over blocks, and only ONE size can be drawn. The size shown is
    # the one with the largest EXCESS OVER THE LINE anywhere on the grid -- not the
    # largest absolute accuracy, which is a different and misleading choice: a 0.5830
    # at the 0.75 column is 0.167 BELOW its bound while a 0.5700 at 0.55 is 0.020
    # ABOVE its own. The excess is the quantity the bound is about, so the size that
    # exceeds by the most is the least flattering one to us. Which size wins gets
    # printed rather than chosen silently.
    sizes = sorted(a34["table"]["a4"], key=int)
    worst = max(sizes, key=lambda s: max(
        _k(a34["table"]["a4"][s][cond], c)["max_balanced"] - c
        for cond in CONDS for c in CEILINGS))
    series[f"A4 — segment LR, max over blocks (size {worst})"] = {
        cond: [(c, (v := _k(a34["table"]["a4"][worst][cond], c))["max_balanced"],
                v["se_at_max"]) for c in CEILINGS]
        for cond in CONDS}

    for name, by_cond in series.items():
        for cond, pts in by_cond.items():
            if len(pts) != len(CEILINGS):
                raise SystemExit(f"{name}/{cond}: {len(pts)} of {len(CEILINGS)} cells")

    # ⚠ `bound_check` is written INSIDE the table (arm2_attack_a3a4.py stores it as
    #   table["bound_check"]), so a top-level lookup silently misses and the corrected
    #   thresholds fall back to literals that would be wrong for any other grid size.
    #   Fail loudly on a shape we do not recognise rather than defaulting.
    bc = (a34.get("table", {}) or {}).get("bound_check") or a34.get("bound_check") or {}
    if not bc:
        print("⚠ this JSON carries no `bound_check` block (it predates it, or was "
              "written by a --from-json re-table, which does not persist one), so the "
              "corrected thresholds below are the values the shipped run printed and "
              "are NOT read from the data.")
    return {"series": series, "a4_size": worst,
            "z_a3": bc.get("z_crit_a3", 2.82), "z_a4": bc.get("z_crit_a4", 4.51),
            "n_per_arm": a1.get("trials_per_cell"), "delta": a1.get("delta")}


def _gate_bars_reach_below(series: dict, z_a3: float, z_a4: float) -> list:
    """⚠⚠ THE GATE THIS FIGURE EXISTS UNDER (Dylan, 2026-10-29).

    A point above the diagonal whose 95% bar does NOT reach below it would be drawn as
    an unqualified exceedance. That must stop the run, not appear in a PDF. Degenerate
    cells (se = 0) are exempt and reported separately -- they carry no bar because they
    carry no information.
    """
    above, bad, degen = [], [], []
    for name, by_cond in series.items():
        z_crit = z_a4 if name.startswith("A4") else z_a3
        for cond, pts in by_cond.items():
            for c, bal, se in pts:
                verdict, z, lo, _hi = classify_excess(bal, se, c, z_crit)
                tag = f"{name} / {cond} @ {c:.2f}: {bal:.4f} vs {c:.4f}"
                if verdict == "within":
                    continue
                if verdict == "degenerate":
                    degen.append(tag)
                elif verdict in ("viol", "bar_clears"):
                    bad.append(f"{tag}, 95% lo {lo:.4f} > line, z {z:+.2f} [{verdict}]")
                else:
                    above.append(f"{tag}, 95% lo {lo:.4f}, z {z:+.2f}")
    if bad:
        raise SystemExit(
            "⚠⚠ a point sits above the acc line with its 95% bar NOT reaching below "
            "it. That cannot be drawn without qualification:\n  " + "\n  ".join(bad))
    return above


def fig_accuracy_vs_ceiling(d: dict, path: Path) -> dict:
    series = d["series"]
    fig, axes = plt.subplots(1, len(series), figsize=(14.6, 5.3), sharey=True)
    lim = (0.46, 0.84)
    for ax, (name, by_cond) in zip(axes, series.items()):
        # the permitted region, y <= x, and the line itself
        ax.fill_between(lim, lim, lim[1], color="#d62828", alpha=0.07, lw=0,
                        zorder=0)
        ax.plot(lim, lim, ls="--", lw=1.3, color="#d62828", zorder=1,
                label="`acc_ceiling` (the bound)")
        # ⚠ a small horizontal DODGE only. The three conditions share each release
        # and sit at identical x, so without it their bars overplot and a reader
        # cannot tell which series clears the line. The offset is cosmetic and the
        # tick labels remain the true `acc_ceiling`.
        for off, (cond, pts) in zip((-0.007, 0.0, 0.007), by_cond.items()):
            x = [c + off for c, _b, _s in pts]
            y = [b for _c, b, _s in pts]
            e = [Z95 * s for _c, _b, s in pts]
            degen = [s == 0 for _c, _b, s in pts]
            ax.errorbar(x, y, yerr=e, fmt=MARK[cond], ms=4.6, lw=1.1, capsize=2.6,
                        color=COL[cond], label=cond, zorder=3, alpha=0.95)
            # ⚠ a zero-SE point is a CONSTANT ANSWER, not a precise one. It is ringed
            # so the eye does not read a bar-less marker as a tight measurement.
            if any(degen):
                ax.scatter([xi for xi, g in zip(x, degen) if g],
                           [yi for yi, g in zip(y, degen) if g],
                           s=118, facecolors="none", edgecolors=COL[cond],
                           lw=1.0, zorder=4)
        ax.set_xlim(*lim); ax.set_ylim(*lim)
        ax.set_xticks(list(CEILINGS))
        ax.set_xticklabels([f"{c:.2f}" for c in CEILINGS], fontsize=8)
        ax.set_title(name, fontsize=9.5)
        ax.set_xlabel("`acc_ceiling` (PathGuard's knob) — and the bound", fontsize=8.5)
        ax.grid(alpha=0.18, lw=0.6)
    axes[0].set_ylabel("balanced accuracy, 95% bar", fontsize=9)
    axes[0].legend(fontsize=7.4, loc="upper left", framealpha=0.9)
    fig.suptitle("ARM 2 — measured attacker accuracy against the exact acc line, "
                 "chr21, share 1.0, K = 88", fontsize=11)
    fig.text(0.5, 0.012,
             "Each target runs at ITS OWN eta_tau, solved so its measured acc_ceiling "
             "lands on the column label — so the bound is the DIAGONAL, not a "
             "horizontal. Bars are 1.96 x SE "
             f"(n = {d['n_per_arm']} per arm per cell). Shaded = forbidden. "
             "Pinsker is an inequality, so the true ceiling lies BELOW the dashed line "
             "and the looseness flatters us. Ringed markers have se = 0 — a constant "
             "answer, not a precise one. Conditions share each release and are never "
             "pooled.",
             ha="center", fontsize=7.1, wrap=True)
    fig.tight_layout(rect=(0, 0.115, 1, 0.94))
    fig.savefig(path.with_suffix(".pdf")); fig.savefig(path.with_suffix(".png"), dpi=220)
    plt.close(fig)
    return {"panels": list(series), "a4_block_size": d["a4_size"]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)

    d = load()
    above = _gate_bars_reach_below(d["series"], d["z_a3"], d["z_a4"])
    meta = fig_accuracy_vs_ceiling(d, a.out / "fig_acc_vs_ceiling")

    print("ARM 2 — measured accuracy against the acc line. Metrics are ROWS, "
          "`acc_ceiling` is the COLUMN.\n")
    print(f"⚠ The bound is the DIAGONAL y = x: in ceiling mode the column label IS "
          f"the measured `acc_ceiling`.")
    print(f"⚠ A4's panel shows block size {d['a4_size']} — the size with the largest "
          f"EXCESS OVER THE LINE anywhere on the grid, which is the least flattering "
          f"choice to us (not the largest absolute accuracy, which sits at a high "
          f"ceiling and is far BELOW its own bound). Its bar is the per-block 95%, "
          f"ANTI-CONSERVATIVE for a max; the simultaneous threshold is "
          f"z = {d['z_a4']:.2f}.")
    print(f"\n{len(above)} cell(s) sit ABOVE the line by point estimate, every one "
          f"with a 95% bar reaching below it:")
    for m in above:
        print(f"    {m}")
    # ⚠⚠ DERIVED, never asserted. This sentence read "the 0.50 and 0.55 columns do
    #   not [resolve], which is where every one of those exceedances sits" -- true of
    #   A3/A4 alone, and FALSE of the combined list it was printed under, which
    #   contains an A1 cell at the 0.60 column.
    thr = 0.5 + d["delta"]
    unres = [c for c in CEILINGS if c <= thr]
    at = sorted({float(m.split("@ ")[1].split(":")[0]) for m in above})
    inside = [c for c in at if c <= thr]
    print(f"\n⚠ delta = {d['delta']:.4f}, so a cell resolves only where the ceiling "
          f"exceeds {thr:.4f} — the "
          f"{', '.join(f'{c:.2f}' for c in unres)} column(s) do not.")
    if at:
        print(f"⚠ The points above the line sit at the "
              f"{', '.join(f'{c:.2f}' for c in at)} column(s); "
              + (f"all of them are inside the unresolved range."
                 if len(inside) == len(at) else
                 f"{len(at) - len(inside)} of {len(at)} are at columns the experiment "
                 f"DOES resolve ("
                 f"{', '.join(f'{c:.2f}' for c in at if c > thr)}), so for those the "
                 f"measurement is the statement and not the ceiling."))
    meta["above_the_line"] = above
    (a.out / "figures.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(f"\nfigure -> {a.out / 'fig_acc_vs_ceiling'}.pdf/.png")
    print(f"JSON   -> {a.out / 'figures.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
