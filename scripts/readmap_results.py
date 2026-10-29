#!/usr/bin/env python
"""READ MAPPING ON THE RELEASED GRAPH — PanMixer's Fig 5A quantity, on our releases.

WHAT THIS MEASURES, and it is UTILITY, not privacy. A pangenome graph exists to be a
reference other people align reads against. If substituting an obfuscated genome into
the released graph degrades alignment for an unrelated donor, the graph is worse at its
job. So: map five EXTERNAL donors' chr21 reads against a graph built from our release
and report what fraction of aligned reads come back perfect, gapless, and at maximum
mapping confidence.

⚠⚠ THE COMPARISON IS THE DELTA, NOT THE LEVEL. Our absolute level cannot be compared to
their published 77.83 / 95.62 / 77.01: we run vg 1.68 against their 1.36, and that
offset was measured at +1.2 / +2.2 points on 2026-09-18 and is still unexplained after
the reference-build hypothesis was tested and REFUTED. What IS comparable is how far
each tool's obfuscated graph falls below its OWN baseline. Theirs falls 0.01 points.

⚠ WHY THE BASELINE WAS RE-RUN rather than reused from 2026-09-18: every graph here,
baseline included, is built through ONE code path -- the same npy -> VCF conversion,
the same bcftools chain, the same `vg autoindex`. That matters more than it looks: the
converter's round trip is slightly lossy (280,989 records against the 287,147 the
September run got building straight from `pangenome.vcf.gz`), so reusing the old
baseline would have put a ~6,000-record confound inside every delta.

⚠ DONORS ARE THE PAPER'S FIVE (HG00138 EUR, HG00635 EAS, HG01112 AMR, HG02698 SAS,
NA18853 AFR), per its Methods. Their `constants.py:42` runs HG01600 in place of
NA18853, and NA18853 (AFR) is the LOWEST donor on every metric, so matching their code
would move our absolute numbers UP and FURTHER from their published figures. Deltas are
unaffected by the choice. See `docs/data.md` under `external/PanMixer/read_fastqs`.

⚠ METRIC DEFINITIONS ARE THEIRS, verbatim from
`external/PanMixer/plot_util/plot_helpers.py:15-19` -- each is a fraction of
`total_aligned`, NOT of reads submitted. A mapping RATE proper (aligned / submitted) is
reported by NEITHER tool.

    python scripts/readmap_results.py
"""

from __future__ import annotations

import argparse
import json
import re
import statistics as st
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DONORS = ("HG00138", "HG00635", "HG01112", "HG02698", "NA18853")
# ⚠⚠ TWO GRIDS, AND THEY ARE NOT POOLABLE. `eta` is the 2026-10-26 run, which was
#    CROSS-BUILD on PanMixer's stale 340,824 tree; `acc` is the reported 0.50-0.80
#    ceiling grid, run one-build on our 340,849 axis by
#    `scripts/readmap_acc_run.sh`. Selected by --grid; `eta` stays the default so the
#    registered 2026-10-26 table re-emits byte-identically.
GRIDS = {
    "eta": {
        "dir": REPO / "logs" / "B_readmap_eta",
        "knob": "eta_tau",
        # "base" is the unobfuscated graph; the rest are HG00438's release at each eta.
        "graphs": ("base", "eta_0", "eta_250", "eta_500", "eta_1262", "eta_2821",
                   "eta_5642"),
        "lbl": {"base": "baseline (unobfuscated)", "eta_0": "0", "eta_250": "250",
                "eta_500": "500", "eta_1262": "1262", "eta_2821": "2821",
                "eta_5642": "5642"},
        "build": "CROSS-BUILD (PanMixer's 340,824-site tree, marked STALE)",
    },
    "acc": {
        "dir": REPO / "logs" / "B_readmap_acc",
        "knob": "acc_ceiling",
        "graphs": ("base", "acc_0p50", "acc_0p55", "acc_0p60", "acc_0p65", "acc_0p70",
                   "acc_0p75", "acc_0p80"),
        # ⚠ 0.50 is `eta` = 0 -- a SECOND control, not a privacy point.
        "lbl": {"base": "baseline (unobfuscated)", "acc_0p50": "0.50 [eta = 0]",
                "acc_0p55": "0.55", "acc_0p60": "0.60", "acc_0p65": "0.65",
                "acc_0p70": "0.70", "acc_0p75": "0.75", "acc_0p80": "0.80"},
        "build": "ONE BUILD (our 340,849-site axis)",
    },
}
D = GRIDS["eta"]["dir"]
GRAPHS = GRIDS["eta"]["graphs"]
LBL = GRIDS["eta"]["lbl"]
KNOB = GRIDS["eta"]["knob"]
BUILD = GRIDS["eta"]["build"]
# PanMixer's published pair, for the DELTA comparison only.
THEIR_BASE = (77.83, 95.62, 77.01)
THEIR_OBF = (77.82, 95.61, 77.01)


def parse(f: Path) -> dict:
    """`vg stats -a` output -> the three published fractions.

    ⚠ Fails loudly on a field it cannot find. A silently-missing numerator would make
    a graph look better than it is, which is the one direction that must not pass.
    """
    txt = f.read_text()

    def one(pat: str, what: str) -> float:
        # ⚠ re.M: every pattern here is anchored at a LINE start. Without it `^` only
        # matches the start of the file and all four lookups fail -- which the gate
        # below catches, but as a format error rather than the bug it is.
        m = re.search(pat, txt, re.M)
        if not m:
            raise SystemExit(f"{f.name}: no {what} -- `vg stats -a` format changed?")
        return float(m.group(1))

    aligned = one(r"^Total aligned: (\d+)", "Total aligned")
    if aligned <= 0:
        raise SystemExit(f"{f.name}: Total aligned is {aligned}")
    perfect = one(r"^Total perfect: (\d+)", "Total perfect")
    gapless = one(r"^Total gapless \(softclips allowed\): (\d+)", "Total gapless")
    # the "(N reads)" that follows "max 60" on the Mapping quality line
    mq = one(r"^Mapping quality:.*max 60 \((\d+) reads\)", "MAPQ-60 read count")
    return {"aligned": aligned,
            "perfect": 100.0 * perfect / aligned,
            "gapless": 100.0 * gapless / aligned,
            "mapq60": 100.0 * mq / aligned}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--grid", choices=tuple(GRIDS), default="eta",
                    help="which run to table. `eta` is the 2026-10-26 CROSS-BUILD run "
                         "and stays the default so its registered table re-emits "
                         "unchanged; `acc` is the one-build 0.50-0.80 ceiling grid. "
                         "⚠ The two are NOT poolable and must not share a table.")
    a = ap.parse_args(argv)

    global D, GRAPHS, LBL, KNOB, BUILD
    g_ = GRIDS[a.grid]
    D, GRAPHS, LBL = g_["dir"], g_["graphs"], g_["lbl"]
    KNOB, BUILD = g_["knob"], g_["build"]

    per: dict = {}
    for g in GRAPHS:
        per[g] = {}
        for s in DONORS:
            f = D / f"stats_{g}_{s}.txt"
            if not f.exists():
                raise SystemExit(f"missing {f} -- the matrix is incomplete; refusing "
                                 f"to table a partial fan-out")
            per[g][s] = parse(f)
    n = sum(len(v) for v in per.values())
    if n != len(GRAPHS) * len(DONORS):
        raise SystemExit(f"parsed {n} cells, expected {len(GRAPHS) * len(DONORS)}")

    med = {g: {k: st.median(per[g][s][k] for s in DONORS)
               for k in ("perfect", "gapless", "mapq60")} for g in GRAPHS}

    print("READ MAPPING ON THE RELEASED GRAPH. ARM 1, chr21, target HG00438, "
          f"5 external donors, median over donors. Metrics are ROWS, `{KNOB}` is the "
          f"COLUMN. vg 1.68.0, single-end, `vg autoindex --workflow giraffe`.")
    print(f"⚠ {BUILD}.")
    print("⚠ Each metric is a percentage of ALIGNED reads, per "
          "`external/PanMixer/plot_util/plot_helpers.py:15-19` — not of reads "
          "submitted. ⚠ UTILITY, not privacy. ⚠ Levels are NOT comparable to their "
          "published figures (vg 1.68 vs 1.36); the DELTA is.\n")

    hdr = ["metric", "baseline"] + [LBL[g] for g in GRAPHS[1:]]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))
    for k, lab in (("perfect", "`perfect` %"), ("gapless", "`gapless` %"),
                   ("mapq60", "`mapq60` %")):
        print(f"| {lab} | " + " | ".join(
            f"{med[g][k]:.3f}" for g in GRAPHS) + " |")
    print("| reads aligned (median) | " + " | ".join(
        f"{st.median(per[g][s]['aligned'] for s in DONORS):,.0f}" for g in GRAPHS)
        + " |")

    print("\n**THE DELTA AGAINST OUR OWN BASELINE** — the only comparable quantity.\n")
    hdr2 = ["metric"] + [LBL[g] for g in GRAPHS[1:]] + ["worst", "PanMixer's own delta"]
    print("| " + " | ".join(hdr2) + " |")
    print("|" + "---|" * len(hdr2))
    for i, (k, lab) in enumerate((("perfect", "`perfect` delta, pts"),
                                  ("gapless", "`gapless` delta, pts"),
                                  ("mapq60", "`mapq60` delta, pts"))):
        ds = [med[g][k] - med["base"][k] for g in GRAPHS[1:]]
        worst = min(ds)
        print(f"| {lab} | " + " | ".join(f"{d:+.3f}" for d in ds)
              + f" | **{worst:+.3f}** | {THEIR_OBF[i] - THEIR_BASE[i]:+.2f} |")

    allds = [med[g][k] - med["base"][k]
             for g in GRAPHS[1:] for k in ("perfect", "gapless", "mapq60")]
    print(f"\nworst degradation anywhere on the grid: **{min(allds):+.3f} points**; "
          f"largest movement in either direction: **{max(allds, key=abs):+.3f}**.")
    print("⚠ A positive delta is the obfuscated graph aligning BETTER than the "
          "baseline, which at this magnitude is not a gain.")
    print("⚠ THE ONE MOVEMENT ABOVE NOISE IS GRAPH-LEVEL, NOT DONOR-LEVEL, and the "
          "distinction matters. At `eta_tau` = 500 ALL FIVE donors drop together "
          "(-0.038 to -0.063 pts each), so that graph really does align slightly "
          "worse for everyone — it is not one noisy donor. But it is NOT monotone in "
          "the knob (1262, 2821 and 5642 return to -0.001/-0.002), so it is "
          "draw-to-draw variation in WHICH variants the released row happens to "
          "carry, not an effect of `eta_tau`. One release per point is the limit "
          "here: with n = 1 graph per column this spread cannot be separated from "
          "sampling, and quoting -0.063 as the knob's cost would overstate it.")
    print("⚠ One target (HG00438) substituted into a 44-member graph, so the released "
          "row is 1 of 44: this measures the cost of obfuscating ONE member, exactly "
          "as their Fig 5A does, and must not be read as the cost of obfuscating all.")
    print("⚠ `e_differ`-scale context: the substituted row differs from the original "
          "in ~137,000-139,000 of 681,648 cells at every eta, so the graph really is "
          "being perturbed; a flat result is NOT a null edit.")

    if a.json:
        a.json.write_text(json.dumps(
            {"target": "HG00438", "donors": list(DONORS), "graphs": list(GRAPHS),
             "per_donor": per, "median": med,
             "their_baseline": THEIR_BASE, "their_obfuscated": THEIR_OBF}, indent=1)
            + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
