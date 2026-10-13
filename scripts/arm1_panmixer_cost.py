#!/usr/bin/env python
"""`utility_loss_panmixer` on OUR arm 1 releases, by PanMixer's own formula (decision M).

⚠⚠ DECISION M (Dylan, 2026-10-12): the head-to-head runs PanMixer's evaluator UNMODIFIED,
so this metric uses THEIR shipped support over all **K = 88** -- the target's own two
haplotypes included -- for our released arrays and theirs alike. That puts the same
denominator convention on both sides. Our RELEASE is still fully leave-one-out at K = 86
(decision E); only this metric's scoring support is K = 88. Never put a K = 88-scored
column beside a K = 86 one.

THE FORMULA, transcribed from external/PanMixer/tools/panmixer/obfuscate.py:207-211 and
279-290:

    not_missing[v] = sum over ALL 88 haplotypes of (allele != -1)
    utility_loss[v] = 1 / not_missing[v], or 0 where not_missing[v] == 0
    cost = sum over altered (v, h) cells of utility_loss[v]

⚠⚠ AND THEY DO NOT MASK -1 BEFORE COMPARING -- `new_block != orig_block` is a raw
inequality. That is their convention and we match it, but it does NOT mean the metric is
symmetric across the two tools, and this script reports the decomposition that shows why:

  * THEIR released object is `original_haplotypes.copy()` with SELECTED BLOCKS replaced
    (`obfuscate.py:330-340`), so it PRESERVES -1 everywhere outside those blocks.
  * OURS is a complete mosaic: `ffbs.release` raises if a released path carries -1 at all,
    because RENORMALISE means the chain only ever occupies a donor with a call.

So every cell where the TARGET is uncalled is altered in our output by construction, while
in theirs it is altered only if it falls in a selected block. Those cells are a small share
of cells but a large share of the cost, because an uncalled position is usually a
low-support one and `1/support` is large there. ⚠ **THIS IS NOT A DEFECT IN PANMIXER.** It
is a difference in what the two mechanisms emit -- a path release is complete by
construction -- and the honest response is to report the headline under their exact
convention AND the decomposition beside it, never to silently mask -1 on one side only.

⚠ SCOPE. Our release lives on the CHAIN axis (332,521 of chr21's 340,849 records); the
off-chain positions are not released at all. Their run covers their own build's axis. So
raw costs are NOT comparable across the two and only SHARES of each side's own total are
-- which is what `docs/data.md` has said about this comparison since 2026-10-07.

    python scripts/arm1_panmixer_cost.py --scope chr21 --json logs/B_arm1_pm_cost.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

from cohort import loo                                     # noqa: E402

PINNED_TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
AUTOSOMES = loo.AUTOSOMES
ARM1 = REPO / "data" / "releases" / "arm1"


def panmixer_weights(chrom: str):
    """`utility_loss[v]` exactly as `obfuscate.py:207-211` builds it, at K = 88."""
    hap = np.asarray(np.load(REPO / "data" / "cohort" / chrom / "haplotypes.npy",
                             mmap_mode="r"))
    nm = (hap != -1).sum(axis=0)                 # == np.sum(pangenome != -1, axis=(0,2))
    ul = np.zeros(nm.shape[0], dtype=np.float64)
    ok = nm > 0
    ul[ok] = 1.0 / nm[ok]
    return hap, ul


def released_for(target: str, tau: float, chrom: str) -> Path:
    """The one stored part for this (target, `tau`, chromosome)."""
    d = ARM1 / target / f"tau_{tau:.2f}".replace(".", "p")
    hits = sorted(d.glob(f"*chrom{chrom}_*_{chrom}/released_path.npy"))
    if len(hits) != 1:
        raise AssertionError(
            f"{target}/{tau:g}/{chrom}: {len(hits)} released_path.npy, expected exactly "
            f"1. Two would mean constraint 5's one-release-per-genome was violated.")
    return hits[0]


def cost_one(target: str, taus, chroms) -> dict:
    out = {"target": target, "cells": {}}
    acc = {t: {"cost": 0.0, "cost_called": 0.0, "cost_uncalled": 0.0,
               "total": 0.0, "altered": 0, "altered_uncalled": 0, "cells": 0,
               "uncalled": 0} for t in taus}
    for chrom in chroms:
        hap, ul = panmixer_weights(chrom)
        # ⚠⚠ THE PER-TARGET CHAIN, NOT arm 2's `in_chain.npy`. This read the stored
        # arm 2 mask until 2026-10-13 and the shape assertion below caught it on 24 of
        # 44 targets: arm 1's chain is per-target and is NOT a subset of arm 2's --
        # `untraversable_block` is not monotone in the donor panel, so HG00621's chr21
        # chain is 332,518 against arm 2's 332,521 while HG02630's is LONGER
        # (logs/B_loo_chain_survey.md). Using arm 2's axis here silently diffed the
        # release against the wrong positions, and the targets that happened to PASS
        # were exactly those whose chain equals arm 2's -- a biased subset.
        ic = loo.loo_chain(REPO / "data" / "cohort" / chrom, target,
                           haplotypes=hap).in_chain
        rows = list(loo.sample_index(np.load(
            REPO / "data" / "cohort" / chrom / "haplotype_ids.npy",
            allow_pickle=True))[target])
        orig = hap[rows][:, ic].T                 # (T, 2) the target's TRUE haplotypes
        ul_c = ul[ic]
        uncalled = orig < 0
        for tau in taus:
            rel = np.load(released_for(target, tau, chrom))
            if rel.shape != orig.shape:
                raise AssertionError(
                    f"{target}/{chrom}: released {rel.shape} against target "
                    f"{orig.shape}; these are not the same axis")
            # ⚠ RAW inequality, no -1 masking: their convention (decision M)
            alt = rel != orig
            w = alt * ul_c[:, None]
            a = acc[tau]
            a["cost"] += float(w.sum())
            a["cost_uncalled"] += float(w[uncalled].sum())
            a["cost_called"] += float(w[~uncalled].sum())
            a["total"] += float(2.0 * ul_c.sum())
            a["altered"] += int(alt.sum())
            a["altered_uncalled"] += int(alt[uncalled].sum())
            a["cells"] += int(alt.size)
            a["uncalled"] += int(uncalled.sum())
        del hap, ul, orig
    for tau in taus:
        a = acc[tau]
        tot = a["total"] or float("nan")
        out["cells"][str(tau)] = {
            **a,
            "share": a["cost"] / tot,
            "share_called_only": a["cost_called"] / tot,
            "uncalled_cost_fraction": a["cost_uncalled"] / a["cost"] if a["cost"] else float("nan"),
            "altered_share": a["altered"] / a["cells"] if a["cells"] else float("nan"),
        }
    return out


def _one(args):
    try:
        return cost_one(*args)
    except Exception as e:
        return {"target": args[0], "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1200:]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scope", choices=("chr21", "genome"), default="chr21",
                    help="chr21 MATCHES PanMixer's run; genome is ours in full")
    ap.add_argument("--targets", nargs="+", default=None)
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--from-json", type=Path, default=None)
    a = ap.parse_args(argv)

    if a.from_json:
        d = json.loads(a.from_json.read_text())
        print(f"re-reporting {len(d['results'])} target(s) from {a.from_json}. "
              f"NOTHING recomputed.")
        _report(d["results"], [float(t) for t in d["taus"]], d["scope"])
        return 0

    chroms = ["chr21"] if a.scope == "chr21" else list(AUTOSOMES)
    members = sorted(loo.sample_index(np.load(
        REPO / "data" / "cohort" / "chr1" / "haplotype_ids.npy", allow_pickle=True)))
    targets = a.targets or members

    print(f"`utility_loss_panmixer` on ARM 1's releases, by PanMixer's formula at "
          f"K = 88 (decision M).")
    print(f"scope {a.scope} ({len(chroms)} chromosome(s)), {len(targets)} target(s) x "
          f"{len(a.taus)} `tau` = {len(targets) * len(a.taus)} cells.")
    print(f"⚠ chr21 scope MATCHES PanMixer's run, which was HG00438 on chr21 from their "
          f"OWN build; raw costs are not comparable across builds, only shares.\n")

    t0, out = time.time(), []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_one, (t, a.taus, chroms)): t for t in targets}
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            out.append(r)
            tag = "❌ " + r["error"] if "error" in r else "✅"
            print(f"  [{n}/{len(targets)}] {tag} {r['target']} "
                  f"({time.time() - t0:.0f}s)", flush=True)

    errs = [r for r in out if "error" in r]
    ok = [r for r in out if "error" not in r]
    assert len(ok) + len(errs) == len(targets), "results lost between submit and collect"
    if ok:
        _report(ok, a.taus, a.scope)
    if a.json:
        a.json.write_text(json.dumps(
            {"arm": "arm1", "scope": a.scope, "support_k": 88, "taus": a.taus,
             "results": out}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    if errs:
        for r in errs[:3]:
            print(f"\n❌ {r['target']}\n{r['traceback']}")
        raise SystemExit(f"INCOMPLETE: {len(errs)} target(s) failed.")
    return 0


def _report(rows, taus, scope) -> None:
    assert rows, "no targets measured; the report would be silent"
    n = len(rows)
    print(f"\n`utility_loss_panmixer` on ARM 1, scope {scope}, support K = 88 "
          f"(decision M), n = {n} held-out individual(s). PanMixer's own formula.\n")
    hdr = ["metric"] + [f"`tau` = {t:g}" for t in taus]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))

    def row(label, key, fmt=".6f", agg=np.median):
        vals = [[r["cells"][str(t)][key] for r in rows] for t in taus]
        print(f"| {label} | " + " | ".join(f"{agg(v):{fmt}}" for v in vals) + " |")

    # ⚠ ONE DRAW, so this is a SAMPLED quantity like `utility_retained`, not an exact
    # expectation. Its variation across `tau` is DRAW NOISE, not a response to `tau`:
    # each cell is an independent draw with its own stream. The exact counterpart is
    # computable -- the cost is a linear functional of the marginals, with weights
    # 1/support(v) in place of `beta_v` -- but it is not computed here.
    row("`utility_loss_panmixer` * [median, 1 draw]", "cost", ",.4f")
    row("its share of our own total [median]", "share", ".6f")
    # ⚠⚠ NOT `target_fidelity`, and the names must not be conflated.
    # `docs/terms.md:target_fidelity` is the differ rate over the target's CALLED cells;
    # this row's denominator is ALL cells and its numerator includes target-uncalled ones.
    # On HG00438 at `tau` = 0 they read 0.191513 and 0.160166 -- the same mechanism, two
    # denominators. This row exists because PanMixer's own "changed=X/Y" may use the
    # all-cells convention, which would make THIS the comparable one.
    row("altered cells / ALL cells [NOT `target_fidelity`]", "altered_share", ".6f")
    row("cells altered [median]", "altered", ",.0f")
    row("our total possible cost [median]", "total", ",.4f")

    print(f"\n⚠⚠ **THE DECOMPOSITION THAT MAKES THIS INTERPRETABLE.** Their formula does "
          f"not mask -1, and our release is COMPLETE by construction while theirs is the "
          f"original with selected blocks replaced, so every target-uncalled cell is "
          f"altered in ours and only some are in theirs:")
    row("  cost from target-CALLED cells [median]", "cost_called", ",.4f")
    row("  cost from target-UNCALLED cells [median]", "cost_uncalled", ",.4f")
    row("  uncalled share OF THE COST [median]", "uncalled_cost_fraction", ".4%")
    row("  share, called cells only [median]", "share_called_only", ".6f")
    u = np.median([rows[0]["cells"][str(taus[0])]["uncalled"]] if len(rows) == 1 else
                  [r["cells"][str(taus[0])]["uncalled"] for r in rows])
    c = np.median([r["cells"][str(taus[0])]["cells"] for r in rows])
    print(f"\n⚠ Target-uncalled cells are only **{u / c:.4%}** of cells "
          f"({u:,.0f} of {c:,.0f}) but carry the cost share above, because an uncalled "
          f"position is usually a LOW-SUPPORT one and `1/support` is large there.")
    print(f"⚠⚠ **NOT A DEFECT IN PANMIXER.** Its goal is donor protection and we claim no "
          f"defect; this is a difference in what the two mechanisms EMIT. A path release "
          f"is complete by construction (`ffbs.release` refuses a released path carrying "
          f"-1), while theirs preserves the original outside the blocks it swaps "
          f"(`obfuscate.py:330-340`). Report the headline under their exact convention "
          f"AND this decomposition; never mask -1 on one side only.")
    print(f"⚠ Raw costs are NOT comparable across the two builds -- their total was "
          f"9288.85 on their axis against ours above on ours, and our release covers the "
          f"CHAIN only. Shares of each side's own total are the quotable figures.")
    print(f"⚠ ARM 1, K = 86 release, K = 88 SCORING support. Never place a K = 88-scored "
          f"column beside a K = 86 one.")


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
