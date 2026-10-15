#!/usr/bin/env python
"""`donor_coherence`: is a release a MOSAIC of real haplotypes, or a patched copy of the target?

WHY THIS EXISTS. Both tools alter about the same share of cells (~16-17%) and both draw
their replacement alleles from real donors at ~99% of altered positions, so neither of
those measurements separates them. `ld_loss` says they differ by 7.4x, which means the
difference is in the COMBINATIONS rather than the individual alleles. This measures the
combinations directly.

THE MEASUREMENT. In a window of W consecutive positions, take the largest share of that
window explained by any SINGLE donor haplotype. A Li-Stephens mosaic is locally ONE donor,
so its best donor explains nearly the whole window. A copy of the target with scattered
cells flipped is explained by no pool donor, because the target is excluded from the pool.

⚠ THE TARGET'S OWN HAPLOTYPE IS THE REFERENCE, not 100%. A real human genome is itself only
~98.6% explained by the best of 86 other haplotypes over 200 positions -- that gap IS human
haplotype sharing. Reading a release against 100% rather than against the real genome would
make every release look worse than it is.

⚠ THE TARGET IS EXCLUDED FROM THE DONOR POOL (43 others x 2 = 86 haplotypes). Leaving it in
would let a patched copy of the target match itself and score ~100%, which measures nothing.

⚠ NOT A PRIVACY METRIC AND NOT AN ENDPOINT. High coherence means the release looks like a
real haplotype; it says nothing about WHOSE. A release can be perfectly coherent and carry
zero target information -- that is exactly what `tau` = 0 is.

    python scripts/arm1_donor_coherence.py [--window 200] [--windows 2000]
"""

from __future__ import annotations

import argparse
import json
import statistics as st
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
PM = REPO / "external" / "PanMixer" / "starting_data" / "chr21"
SWEEP = Path("/data/ds85/RanPanMixer/runs/capacity_sweep")
ROOT = Path("/data/ds85/RanPanMixer/runs/beagle_attack/root")
TARGET = "HG00438"
CAPS = ["0p01", "0p05", "0p1", "0p2", "0p3", "0p5", "0p7", "0p9"]
TAUS = ["0p00", "0p05", "0p10", "0p25", "0p50", "0p90", "0p99"]


def coherence(hap_2col, pool, starts, W, min_called=0.5):
    """Mean over windows of (best single donor's match share). Both haplotype columns."""
    out = []
    for col in (0, 1):
        h = hap_2col[:, col]
        for s in starts:
            sl = slice(s, s + W)
            hw, P = h[sl], pool[:, sl]
            valid = (hw >= 0) & (P >= 0)
            n = valid.sum(axis=1)
            ok = n >= W * min_called
            if not ok.any():
                continue
            m = ((P == hw[None, :]) & valid).sum(axis=1)
            out.append(float((m[ok] / n[ok]).max()))
    return st.mean(out), len(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--windows", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    orig = np.load(PM / "pangenome.npy")
    subs = np.load(PM / "pangenome_subjects.npy").tolist()
    ti = subs.index(TARGET)
    others = np.delete(orig, ti, axis=0)
    pool = np.concatenate([others[:, :, 0], others[:, :, 1]], axis=0)
    T = pool.shape[1]
    rng = np.random.default_rng(a.seed)
    starts = rng.choice(T - a.window, size=a.windows, replace=False)

    print(f"`donor_coherence`, chr21, target {TARGET}, ARM 1.")
    print(f"{a.windows:,} random windows of {a.window} positions; donor pool {pool.shape[0]} "
          f"haplotypes with the target EXCLUDED.\n")
    ref, n = coherence(orig[ti], pool, starts, a.window)
    print(f"**REFERENCE** -- the target's OWN true haplotypes: **{ref:.4%}** ({n:,} windows)")
    print("  ⚠ This, not 100%, is what a real human genome scores. The shortfall is human "
          "haplotype sharing.\n")

    print("| release | `donor_coherence` | vs the real genome |")
    print("|---|---|---|")
    rows = {r["row_id"]: r for r in json.loads((ROOT / "manifest.json").read_text())["rows"]}
    ours = {}
    for i, r in rows.items():
        if r["subject"] == TARGET and r["tau"] in TAUS:
            h = np.load(ROOT / "experiments/exp_1/data/chr21" / str(i) / "new_haplotypes.npy")
            ours[r["tau"]] = coherence(h, pool, starts, a.window)[0]
    for t in TAUS:
        if t in ours:
            print(f"| ours, `tau` = {t.replace('p', '.')} | {ours[t]:.4%} | {ours[t] - ref:+.4%} |")
    theirs = {}
    for c in CAPS:
        h = np.load(SWEEP / f"cap{c}" / "chr21" / "new_haplotypes.npy")
        theirs[c] = coherence(h, pool, starts, a.window)[0]
        print(f"| theirs, `capacity` = {c.replace('p', '.')} | {theirs[c]:.4%} | "
              f"{theirs[c] - ref:+.4%} |")

    ov = list(ours.values())
    sat = [theirs[c] for c in ("0p3", "0p5", "0p7", "0p9")]
    print(f"\nours span {min(ov):.4%}..{max(ov):.4%}; theirs at saturation "
          f"{min(sat):.4%}..{max(sat):.4%}")
    print(f"⚠ NOT an endpoint and NOT a privacy metric: coherence says the release looks like "
          f"a real haplotype, never WHOSE. A `tau` = 0 release is coherent and carries zero "
          f"target information.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
