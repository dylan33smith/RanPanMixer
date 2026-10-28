#!/usr/bin/env python
"""Our 1-of-2 discrimination attack run on PANMIXER'S OWN RELEASES, same 44 people.

Spec: `docs/memory.md` 2026-10-27 (B-ATK-arm1-spec) — "run A1 and A3 on Arm 1's 308 chr21
releases, and on PanMixer's, since A1 is `lrt_pair` and well-posed in both arms". The
POSING is well posed in both arms; A1's STATISTIC is not available for their mechanism,
and that distinction is the first thing this script reports rather than papers over.

⚠⚠ A1 IS `n/a` HERE, WITH A REASON, NOT OMITTED. A1 is the Neyman-Pearson optimal test
for OUR release distribution: `Lambda = eta[u(p,y) - u(q,y)] - [log Z_p - log Z_q]`,
which exists in closed form only because `Q_p(y) ∝ R_D(y) exp(eta u(p,y))`. PanMixer's
release is the output of a constrained optimiser followed by a resampling step; it has no
`eta`, no `Z`, and no tractable likelihood, so no Neyman-Pearson statistic can be written
for it. Substituting a heuristic and calling it A1 would compare the OPTIMAL attacker on
one tool against an arbitrary one on the other, and the null on their side would then be
a statement about our choice of heuristic. So A1 runs on our releases only, and this
script runs A3 — which is mechanism-agnostic because it is PanMixer's OWN score.

⚠⚠ WHAT IS AND IS NOT COMPARABLE BETWEEN THE TWO TABLES.
  COMPARABLE: the ACCURACIES, TPR and TNR. They are dimensionless and both sides answer
    the identical question on the identical (target, decoy) pairs — the decoy identities
    are READ FROM our own run's JSON rather than re-derived, so the pairs cannot drift.
  NOT COMPARABLE: the raw A3 SCORES. Our Arm 1 axis is the per-target leave-one-out chain
    (~332,521 positions, 86 donor haplotypes); theirs is their full pangenome axis
    (340,849 positions, 88 haplotypes, nobody held out). The bi-allelic site set and the
    HWE rarity weights therefore differ, so the score MAGNITUDES live on two different
    scales. Never put them in one column.

⚠⚠ THERE IS NO `acc_ceiling` ON THIS SIDE, and that is the asymmetry worth reporting.
Theorem 1 gives our releases a per-pair bound that a measured accuracy can be tested
against. PanMixer's privacy notion is its optimiser's gap-score constraint, which bounds
nothing of this form, so a measured accuracy here has nothing to exceed. The row is
printed as `n/a` with that reason. ⚠ This cuts BOTH ways and must be reported as such: a
low measured accuracy on their side is not evidence of a guarantee, and our bound is not
evidence that our measured accuracy is lower.

⚠ D1_extremal IS `n/a` HERE TOO. It is a CONSTRUCTED path that differs from the target at
every called position; it is not a cohort member, so PanMixer never obfuscated it and no
release exists under that hypothesis. Only D2_nearest and D3_random have both arms.

⚠ THE TARGET IS NOT HELD OUT on their side — PanMixer obfuscates a copy of the subject's
own row inside the full 44-member pangenome. So the bi-allelic axis and the HWE weights
are computed from a panel that CONTAINS the target. That is their design, not a defect,
and it is the mirror of our own panel asymmetry (we hold the target out and the decoys
stay in). Both asymmetries are stated with their tables.

    python scripts/pm_attack.py --ours logs/B_arm1_attack/attack.json
"""

from __future__ import annotations

import argparse
import json
import math
import statistics as stats
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src"), str(REPO / "scripts")]

import numpy as np                                             # noqa: E402

import attacks as A                                            # noqa: E402
import arm2_attack_a3a4 as A3M                                 # noqa: E402

SD = Path("/home/ds85/pm_ourbuild/starting_data/chr21")
OBF = Path("/home/ds85/pm_ourbuild/obf_out")                   # seed 123
DRAWS = Path("/home/ds85/pm_ourbuild/obf_draws")               # seeds 124..127
CAPS = ("0p01", "0p05", "0p1", "0p2", "0p3", "0p5", "0p7", "0p9")
CONDS = ("D2_nearest", "D3_random")
DEGENERATE_FLOOR = 0.05


def load_axis() -> tuple:
    """PanMixer's own pangenome, as haplotype rows, plus the A3 site axis on it.

    ⚠ `pangenome.npy` is (n_samples, T, 2); `a3_sites` wants (K, T) haplotype rows. The
    transpose-then-reshape must keep each sample's two strands ADJACENT and in order, or
    the derived `alt_af` would be a permutation of the truth and every weight would be
    wrong while every number stayed plausible. Asserted, not assumed.
    """
    pan = np.load(SD / "pangenome.npy")
    subs = [str(s) for s in np.load(SD / "pangenome_subjects.npy")]
    hap = pan.transpose(0, 2, 1).reshape(pan.shape[0] * 2, pan.shape[1])
    for i in (0, pan.shape[0] // 2, pan.shape[0] - 1):
        assert (hap[2 * i] == pan[i, :, 0]).all(), f"strand 0 misplaced at sample {i}"
        assert (hap[2 * i + 1] == pan[i, :, 1]).all(), f"strand 1 misplaced at {i}"
    A3M._A3.clear()
    s = A3M.a3_sites({"hap": hap})
    assert len(s["mask"]) == pan.shape[1]
    return pan, subs, s


def releases_for(subject: str, cap: str) -> list:
    """All five of PanMixer's draws for one (subject, capacity).

    ⚠ The four extra seeds live under `obf_draws/<s>/cap*/s12[4-7]/`, where the VCFs were
    deleted to save 53 GB but `new_haplotypes.npy` SURVIVED — measured, 1,408 of them. So
    their side has 5 draws per cell, not 1, and the count is ASSERTED per cell rather
    than inferred from the directory listing.
    """
    out = [OBF / subject / f"cap{cap}" / "chr21" / "new_haplotypes.npy"]
    out += [DRAWS / subject / f"cap{cap}" / f"s{s}" / "chr21" / "new_haplotypes.npy"
            for s in (124, 125, 126, 127)]
    missing = [str(p) for p in out if not p.exists()]
    if missing:
        raise SystemExit(f"{subject}/cap{cap}: {len(missing)} draw(s) missing: "
                         f"{missing[:2]}")
    return [np.load(p) for p in out]


def _acc(lp: list, lq: list) -> tuple:
    a, b = np.asarray(lp, float), np.asarray(lq, float)
    tpr = float(((a > 0).sum() + 0.5 * (a == 0).sum()) / len(a))
    tnr = float(((b < 0).sum() + 0.5 * (b == 0).sum()) / len(b))
    bal = A.empirical_accuracy(a, b)
    assert abs(bal - (0.5 * tpr + 0.5 * tnr)) < 1e-12, "accuracy definitions drifted"
    se = 0.5 * math.sqrt(tpr * (1 - tpr) / len(a) + tnr * (1 - tnr) / len(b))
    return bal, se, tpr, tnr


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", type=Path,
                    default=REPO / "logs" / "B_arm1_attack" / "attack.json",
                    help="our own run's JSON; the decoy identities are read from it")
    a = ap.parse_args(argv)
    if not a.ours.exists():
        raise SystemExit(
            f"{a.ours} missing. The decoy identities are READ from our own run rather "
            f"than re-derived, so the two tools are scored on identical (target, decoy) "
            f"pairs; without it this script would pick its own decoys and the "
            f"comparison would silently be between two different pairings.")
    ours = json.loads(a.ours.read_text())
    label = ours["labels"][-1]
    # ⚠ the decoy is a property of the TARGET, not of the acc_ceiling label -- the same
    # pair is used at every label -- so reading it from any one label is sound. Asserted.
    decoy: dict = {}
    for r in ours["results"]:
        for cond in CONDS:
            ids = {r["cells"][lab]["conditions"][cond]["decoy"]
                   for lab in ours["labels"]}
            assert len(ids) == 1, f"{r['target']}/{cond}: decoy moves across labels: {ids}"
            decoy[(r["target"], cond)] = ids.pop()

    pan, subs, sites = load_axis()
    targets = sorted({r["target"] for r in ours["results"]})
    idx = {s: i for i, s in enumerate(subs)}
    miss = [t for t in targets if t not in idx]
    if miss:
        raise SystemExit(f"not in PanMixer's pangenome_subjects.npy: {miss}")

    res: dict = {}
    for cap in CAPS:
        for cond in CONDS:
            lp, lq = [], []
            for t in targets:
                q = decoy[(t, cond)]
                tp, tq = pan[idx[t]], pan[idx[q]]
                for y in releases_for(t, cap):
                    lp.append(A3M.a3_score(y, tp, sites) - A3M.a3_score(y, tq, sites))
                for y in releases_for(q, cap):
                    lq.append(A3M.a3_score(y, tp, sites) - A3M.a3_score(y, tq, sites))
            assert len(lp) == len(lq) == 5 * len(targets), (
                f"cap{cap}/{cond}: {len(lp)}/{len(lq)} draws, expected "
                f"{5 * len(targets)} each")
            res[(cap, cond)] = _acc(lp, lq)
        print(f"  cap{cap} done", file=sys.stderr)

    n = 5 * len(targets)
    print("# A3 ON PANMIXER'S OWN RELEASES — the same 44 people, the same decoy pairs\n")
    print(f"chr21. Their unmodified mechanism on OUR build, {len(targets)} subjects x "
          f"{len(CAPS)} capacities x 5 draws per hypothesis; n = {n} per arm per cell.")
    print(f"Decoy identities READ FROM `{a.ours}` so both tools are scored on identical "
          f"(target, decoy) pairs. A3's site axis is PanMixer's own "
          f"({sites['n_bi']:,} bi-allelic of {pan.shape[1]:,}), 88 haplotypes, "
          f"NOBODY held out.\n")
    print("⚠⚠ **A1 IS `n/a` ON THIS SIDE.** It is the Neyman-Pearson optimum for OUR "
          "release distribution and needs `eta`, `u` and `log Z`; PanMixer's release is "
          "an optimiser output with no tractable likelihood, so no such statistic "
          "exists for it. A heuristic called A1 would compare the optimal attacker on "
          "one tool against an arbitrary one on the other.")
    print("⚠⚠ **`acc_ceiling` IS `n/a` ON THIS SIDE.** Theorem 1 gives our releases a "
          "per-pair bound a measured accuracy can be tested against; PanMixer's privacy "
          "notion bounds nothing of this form. ⚠ That cuts BOTH ways: a low accuracy "
          "here is not evidence of a guarantee, and our having a bound is not evidence "
          "that our measured accuracy is lower.")
    print("⚠ **D1_extremal is `n/a`**: it is a constructed path, not a cohort member, so "
          "PanMixer never obfuscated it and the second hypothesis has no draws.")
    print("⚠ The target is NOT held out on their side — the axis and the HWE weights "
          "come from a panel CONTAINING it. The mirror of our own asymmetry.")
    print("⚠ Raw A3 SCORES are NOT comparable across the two tables (different site "
          "sets and weights). The accuracies are.\n")

    for cond in CONDS:
        print(f"## {cond}\n")
        print("| metric | " + " | ".join(c.replace("0p", "0.") for c in CAPS) + " |")
        print("|---|" + "---|" * len(CAPS))
        rows = {"balanced accuracy*": [], "SE": [], "TPR": [], "TNR": [],
                "`acc_ceiling`*": [], "verdict": []}
        for cap in CAPS:
            bal, se, tpr, tnr = res[(cap, cond)]
            rows["balanced accuracy*"].append(f"{bal:.6f}")
            rows["SE"].append(f"{se:.6f}")
            rows["TPR"].append(f"{tpr:.4f}")
            rows["TNR"].append(f"{tnr:.4f}")
            rows["`acc_ceiling`*"].append("n/a")
            rows["verdict"].append("DEGENERATE" if min(tpr, tnr) < DEGENERATE_FLOOR
                                   else "no bound to exceed")
        for k, v in rows.items():
            print(f"| {k} | " + " | ".join(v) + " |")
        print()

    print("\n## SYNTHESIS\n")
    for cond in CONDS:
        vals = [res[(c, cond)][0] for c in CAPS]
        lo, hi = min(vals), max(vals)
        sat = res[("0p9", cond)]
        print(f"- **{cond}**: A3's balanced accuracy spans {lo:.6f}..{hi:.6f} across "
              f"their capacity range; at saturation (`capacity` 0.9) it is "
              f"{sat[0]:.6f} +/- {sat[1]:.6f}, TPR {sat[2]:.4f} / TNR {sat[3]:.4f}.")
    print(f"- ⚠ Chance is 0.5 and n = {n} per arm, so SE at chance is "
          f"{0.5 * math.sqrt(2 * 0.25 / n):.6f}: an accuracy inside "
          f"{0.5 - 1.96 * 0.5 * math.sqrt(0.5 / n):.4f}.."
          f"{0.5 + 1.96 * 0.5 * math.sqrt(0.5 / n):.4f} is indistinguishable from "
          f"chance at this n.")
    print("- ⚠ An attacker that fails is not evidence the bound holds (standing "
          "constraint 4), and A3 is NOT optimal for either mechanism. A null here is a "
          "statement about A3, not about PanMixer.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
