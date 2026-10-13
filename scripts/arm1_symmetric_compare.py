#!/usr/bin/env python
"""BOTH releases scored by BOTH metrics, on the SAME positions. The symmetric table.

⚠⚠ WHY THIS EXISTS. `scripts/arm1_headtohead.py` scores each tool on its own axis, so
only SHARES are comparable and our own utility function appears on our side alone. This
restricts to the positions the two builds SHARE -- 339,744 of chr21, 99.68% of each --
and evaluates every metric on every release. One axis, one set of weights, two objects.

⚠ WHAT BECOMES COMPARABLE, AND WHY. `utility_loss_panmixer` and `target_fidelity` were
already comparable. `e_q_u` was not, because PanMixer has no `beta_v` -- it has no utility
function at all, only a cost. But OUR utility is just a weighted match count, and a
weighted match count can be evaluated on ANY released object. So we score THEIR release
with OUR `beta_v` and ours with theirs. **That is not PanMixer reporting a utility; it is
us measuring their output on our axis**, and it must be labelled that way.

⚠⚠ THE WEIGHTS ARE OURS AND THEY ARE INHERITED FROM THEM. `prereg_2026-10-10.md` records
that `w_v = 1/support(v)` "came from PanMixer, where it was a cost weight doing a
different job, and it upweights poorly-called rather than informative positions" -- it is
NOT a rarity weight. On this cohort it is nearly uniform: 98.6% of chain positions sit at
support >= 80 and hold 96.8% of the weight.

⚠ SCOPE. chr21, the shared sub-axis, so neither side's number here equals its number on
its own full axis. Our `beta_v` is a PARTIAL SUM over the shared positions and no longer
sums to chr21's `beta_share`; it is a fixed weighting applied to two objects, which is
what makes the comparison valid, not a budget here.

⚠ ARM 1, K = 86 release, K = 88 scoring support (decision M).

    python scripts/arm1_symmetric_compare.py --json logs/B_arm1_symmetric.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

from cohort import loo                                     # noqa: E402
from mechanism import utility                              # noqa: E402

PINNED_TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
PM_SWEEP = Path("/data/ds85/RanPanMixer/runs/capacity_sweep")
PM_DATA = REPO / "external" / "PanMixer" / "starting_data" / "chr21"
TARGET = "HG00438"


def score(orig, final, w_cost, beta):
    """Every metric, on one (n, 2) original/final pair. ONE function for both tools."""
    alt = final != orig                                   # raw, obfuscate.py:290
    valid = (orig != -1) & (final != -1)                  # obfuscate.py:474-481
    match = (final == orig) & valid
    nv = int(valid.sum())
    return {
        "cost_raw": float((alt * w_cost[:, None]).sum()),
        "cost_masked": float(((alt & valid) * w_cost[:, None]).sum()),
        "changed": int((alt & valid).sum()),
        "valid_cells": nv,
        "target_fidelity": float((alt & valid).sum()) / nv if nv else float("nan"),
        # OUR utility function, evaluated on whichever object is passed in.
        # ⚠ 0.5 per haplotype: diploid utility is the MEAN of the two.
        "our_utility": float(0.5 * (match * beta[:, None]).sum()),
    }


def pm_gt(vcf: Path, sample: str, n: int) -> np.ndarray:
    out = subprocess.run(["bcftools", "query", "-s", sample, "-f", "[%GT]\\n", str(vcf)],
                         capture_output=True, text=True, check=True)
    g = out.stdout.split()
    assert len(g) == n, f"{vcf}: {len(g)} genotypes against {n} sites"
    a = np.full((n, 2), -1, dtype=np.int16)
    for i, s in enumerate(g):
        parts = s.replace("/", "|").split("|")
        for j in (0, 1):
            x = parts[j] if j < len(parts) else "."
            if x not in (".", ""):
                a[i, j] = int(x)
    return a


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    # ---- the shared axis
    d = REPO / "data" / "cohort" / "chr21"
    our_pos = np.load(d / "positions.npy")
    their_pos = np.load(PM_DATA / "pangenome_positions.npy")
    hap = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
    ic = loo.loo_chain(d, TARGET, haplotypes=hap).in_chain     # per-target chain

    shared = np.intersect1d(our_pos[ic], their_pos)
    oi = np.searchsorted(our_pos, shared)                      # index into our site axis
    ti = np.searchsorted(their_pos, shared)                    # index into theirs
    assert np.array_equal(our_pos[oi], shared) and np.array_equal(their_pos[ti], shared)
    # map shared -> position within OUR CHAIN (released_path is chain-indexed)
    chain_of_site = np.full(our_pos.shape[0], -1, dtype=np.int64)
    chain_of_site[ic] = np.arange(int(ic.sum()))
    ci = chain_of_site[oi]
    assert (ci >= 0).all(), "a shared position is outside our chain"
    print(f"shared axis: {len(shared):,} positions of chr21 "
          f"({len(shared)/int(ic.sum()):.4%} of our chain, "
          f"{len(shared)/len(their_pos):.4%} of their axis)")

    # ---- the two originals must be the SAME GENOME on the shared axis, or nothing
    #      downstream is a comparison of mechanisms.
    pg = np.load(PM_DATA / "pangenome.npy", mmap_mode="r")
    subj = list(np.load(PM_DATA / "pangenome_subjects.npy", allow_pickle=True))
    their_orig = np.asarray(pg[subj.index(TARGET)])[ti]
    rows = list(loo.sample_index(np.load(d / "haplotype_ids.npy", allow_pickle=True))[TARGET])
    our_orig = hap[rows][:, oi].T
    both = (our_orig != -1) & (their_orig != -1)
    agree = int((our_orig[both] == their_orig[both]).sum())
    print(f"the two builds' HG00438 agree at {agree:,} of {int(both.sum()):,} "
          f"cells called by both ({agree/int(both.sum()):.6%})")
    disagree = 1.0 - agree / int(both.sum())
    print(f"⚠⚠ THE TWO BUILDS DISAGREE ABOUT THE TARGET ITSELF at "
          f"{int(both.sum()) - agree:,} cells ({disagree:.4%}). Each tool is scored "
          f"against ITS OWN original, so each number is internally sound -- but this "
          f"sets a NOISE FLOOR on any cross-tool comparison, and it is larger than the "
          f"fidelity difference between the two tools.")
    if disagree > 0.05:
        raise SystemExit("the two originals disagree too much to compare mechanisms")

    # ---- weights, on the shared axis
    nm = (hap != -1).sum(axis=0)[oi]                 # K = 88, decision M
    w_cost = np.zeros(len(shared)); w_cost[nm > 0] = 1.0 / nm[nm > 0]
    lc = loo.loo_chain(d, TARGET, haplotypes=hap)
    sh = loo.loo_shares({c: loo.loo_chain(REPO / "data" / "cohort" / c, TARGET).t_chain
                         for c in loo.AUTOSOMES})
    beta_full = utility.beta_v(lc.support, in_chain=ic, share=sh["chr21"],
                               chrom="chr21").beta
    beta = beta_full[ci]
    print(f"our beta_v over the shared axis sums to {beta.sum():.9f} "
          f"(a PARTIAL SUM of chr21's share {sh['chr21']:.9f}); "
          f"their cost weights sum to {2*w_cost.sum():.4f}\n")

    ours, theirs = [], []
    for tau in a.taus:
        td = REPO / "data" / "releases" / "arm1" / TARGET / f"tau_{tau:.2f}".replace(".", "p")
        hits = sorted(td.glob("*chromchr21_*_chr21/released_path.npy"))
        assert len(hits) == 1, f"{tau}: {len(hits)} parts"
        ours.append({"tau": float(tau),
                     **score(our_orig, np.load(hits[0])[ci], w_cost, beta)})
    for cd in sorted(PM_SWEEP.glob("cap*")):
        st = json.loads((cd / "chr21" / "stats.json").read_text())
        fin = pm_gt(cd / "chr21" / "obfuscated.vcf.gz", TARGET, len(their_pos))[ti]
        theirs.append({"capacity": float(st["capacity"]),
                       **score(their_orig, fin, w_cost, beta)})

    _report(ours, theirs, len(shared), beta.sum())
    if a.json:
        a.json.write_text(json.dumps({"arm": "arm1", "scope": "chr21_shared_axis",
                                      "target": TARGET, "n_positions": len(shared),
                                      "ours": ours, "theirs": theirs}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


def _report(ours, theirs, n, beta_sum) -> None:
    print(f"BOTH RELEASES, BOTH METRICS, ONE AXIS. chr21 shared sub-axis, {n:,} "
          f"positions, target {TARGET}, ARM 1 (K = 86 release), scoring support K = 88.\n")
    for name, rows, knob in (("OURS", ours, "tau"), ("THEIRS", theirs, "capacity")):
        print(f"**{name}**\n")
        h = ["metric"] + [f"`{knob}` = {r[knob]:g}" for r in rows]
        print("| " + " | ".join(h) + " |"); print("|" + "---|" * len(h))
        for lab, k, f in (
                ("`target_fidelity` * [differ rate, lower is better]", "target_fidelity", ".6f"),
                ("cells changed", "changed", ",.0f"),
                ("`utility_loss_panmixer` * [masked]", "cost_masked", ".4f"),
                ("our utility, `beta_v`-weighted MATCH", "our_utility", ".9f")):
            print(f"| {lab} | " + " | ".join(f"{r[k]:{f}}" for r in rows) + " |")
        print()
    ou = [r["our_utility"] for r in ours]
    tu = [r["our_utility"] for r in theirs]
    print(f"⚠⚠ **OUR UTILITY FUNCTION SCORES THEIR RELEASE TOO, and that is what makes "
          f"the last row comparable.** It is a `beta_v`-weighted count of cells the "
          f"release got RIGHT, and a weighted count can be evaluated on any object. "
          f"**This is NOT PanMixer reporting a utility** — it has no utility function, "
          f"only a cost — it is us measuring their output with our yardstick.")
    print(f"⚠ Ours spans {min(ou):.9f} to {max(ou):.9f} across the `tau` grid; theirs "
          f"spans {min(tu):.9f} to {max(tu):.9f} across capacity. The ceiling on this "
          f"sub-axis is {beta_sum:.9f}, which is what a release matching the target "
          f"everywhere would score.")
    print(f"⚠ `beta_v` here is a PARTIAL SUM over the shared positions, so it does not "
          f"sum to chr21's `beta_share`. That is fine for a comparison — it is one fixed "
          f"weighting applied to two objects — but it is not a budget at this scope and "
          f"no privacy statement follows from it.")
    print(f"⚠ `tau` and `capacity` are NOT comparable parameters and are never paired.")
    print(f"⚠ Our weighting is INHERITED FROM PANMIXER (prereg: it 'came from PanMixer, "
          f"where it was a cost weight doing a different job, and it upweights "
          f"poorly-called rather than informative positions'). It is not a rarity weight.")


if __name__ == "__main__":
    sys.exit(main())
