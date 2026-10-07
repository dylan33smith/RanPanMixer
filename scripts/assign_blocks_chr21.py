"""Map our PLINK LD blocks onto the site axis -> `blocks_idx.npy`. Model-B input.

OFF THE MAIN PATH, deliberately. `docs/plan.md` keeps the block build as a separate
optional script so the model-B fallback stays one command away without blocks
becoming a dependency of model C, which does not use them.

⚠ Blocks are a SEGMENTATION of the output support, so where they come from is a
standing-constraint-1 question. These come from the 30x panel with the 10 test
targets removed (`scripts/build_blocks_chr21.sh`). Positions in no called block
become their own singleton unit, which is why model B's chain is not simply the
block count.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import store  # noqa: E402

COH = REPO / "data" / "cohort" / "chr21"
det = COH / "blocks.blocks.det"
if not det.exists():
    raise SystemExit(f"{det} missing -- run scripts/build_blocks_chr21.sh first")

rows = np.genfromtxt(det, skip_header=1, usecols=(1, 2, 3, 4), dtype=float)
bp1, bp2 = rows[:, 0].astype(np.int64), rows[:, 1].astype(np.int64)
order = np.argsort(bp1)
bp1, bp2 = bp1[order], bp2[order]

pos = np.load(COH / "positions.npy")
in_chain = np.load(COH / "in_chain.npy")

# a position belongs to the called block whose [BP1, BP2] contains it; else -1
i = np.searchsorted(bp1, pos, side="right") - 1
inside = (i >= 0) & (pos <= np.where(i >= 0, bp2[np.clip(i, 0, None)], -1))
blk = np.where(inside, i, -1).astype(np.int64)

# give every orphan its own id so consecutive orphans do NOT merge into one unit
orph = ~inside
blk[orph] = len(bp1) + np.arange(int(orph.sum()))

meta = store.load_prov(COH / "positions.npy")
store.save(COH / "blocks_idx.npy", blk,
           inputs={"positions": meta["sha256"],
                   "blocks_det": store.sha256_file(det)},
           params={"n_called_blocks": int(len(bp1)),
                   "panel": "1000g_30x_phased, 3192 samples (10 test targets removed)",
                   "plink": "v1.9.0-b.8 --blocks no-pheno-req --snps-only just-acgt",
                   "orphan_policy": "own singleton unit, never merged with a neighbour"},
           axis_digest=meta.get("axis_digest"),
           note="MODEL-B chain unit per site-axis row. Off the main path; model C "
                "does not read this.")

ch = blk[in_chain]
n_units = 1 + int((np.diff(ch) != 0).sum())
sizes = np.diff(np.concatenate(([0], np.flatnonzero(np.diff(ch) != 0) + 1, [len(ch)])))
print(f"called blocks                     {len(bp1):,}")
print(f"chain positions                   {int(in_chain.sum()):,}")
print(f"  inside a called block           {int(inside[in_chain].sum()):,} "
      f"({inside[in_chain].mean():.1%})")
print(f"  orphans (own singleton unit)    {int(orph[in_chain].sum()):,} "
      f"({orph[in_chain].mean():.1%})")
print(f"model-B chain units               {n_units:,}  (model C: {int(in_chain.sum()):,})")
print(f"  unit size median {np.median(sizes):.0f} p90 {np.percentile(sizes,90):.0f} "
      f"max {sizes.max()}")
