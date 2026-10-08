"""Whole-chromosome mismatch density per `tau`, binned — the data behind the map.

Covers ALL chain positions, not a window. 540 bins of about 564 variants each,
density computed PER CELL so it averages to the same share quoted elsewhere.

⚠ A 200-variant window was shown first and that was the wrong picture: 0.07% of the
chromosome, chosen for maximum churn, so both its coverage and its local density were
unrepresentative. Single-variant resolution needs about 3 px per variant, which 540 px
cannot give 304,664 of them, so the whole-chromosome view must bin — and the binning
is the honest way to show it, with a window kept only as a labelled detail.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mechanism import calibration as C   # noqa: E402
from mechanism import ffbs as F          # noqa: E402
from mechanism import hmm as H           # noqa: E402
from mechanism import utility as U       # noqa: E402

TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
NB = 540
COH = REPO / "data" / "cohort" / "chr21"

hap_full = np.load(COH / "haplotypes.npy")
blocks = np.load(COH / "blocks_idx.npy")
in_chain = np.load(COH / "in_chain.npy") & ~H.dead_block_positions(
    hap_full, blocks, np.load(COH / "in_chain.npy"))
hap = hap_full[:, in_chain]
util = U.beta_v(np.load(COH / "support.npy"), in_chain=in_chain)
tc = np.load(REPO / "data" / "targets" / "HG00096" / "chr21" / "path.npy")[in_chain]
pos = np.load(COH / "positions.npy")[in_chain]
known = tc >= 0

out = {"t": list(TAUS), "n": int(in_chain.sum()), "nb": NB, "seed": SEED, "m": {}}
print(f"### Whole-chromosome mismatch map — HG00096, chr21, seed {SEED}\n")
print(f"{int(in_chain.sum()):,} chain positions in {NB} bins "
      f"({in_chain.sum()/NB:.0f} variants per bin)\n")
print("| model | mean density | min bin | max bin | max shift caused by `tau` | "
      "bins holding a flip at `tau`=0.99 |")
print("|---|---|---|---|---|---|")
for nm, blk in (("C", None), ("B", blocks)):
    m = H.from_arrays(hap_full, np.load(COH / "genetic_pos.npy"),
                      np.load(COH / "positions.npy"), in_chain=in_chain, blocks=blk)
    e = np.linspace(0, m.n_positions, NB + 1).astype(int)
    kb = np.array([known[e[b]:e[b + 1]].sum() for b in range(NB)])
    dens = np.zeros((len(TAUS), NB))
    changed_bins = np.zeros((len(TAUS), NB), dtype=np.int64)
    base_any = None
    for i, tau in enumerate(TAUS):
        r = F.release(m, hap, tc, util, tau, seed=SEED, target_id="HG00096")
        mism = (r.released_path != tc) & known
        for b in range(NB):
            dens[i, b] = mism[e[b]:e[b + 1]].sum() / kb[b] if kb[b] else 0.0
        anyd = mism.any(axis=1)
        if base_any is None:
            base_any = anyd
        ch = anyd != base_any
        for b in range(NB):
            changed_bins[i, b] = ch[e[b]:e[b + 1]].sum()
    out["m"][nm] = {
        "d": np.round(dens, 4).tolist(),
        "cb": [[[int(b), int(v)] for b, v in enumerate(row) if v] for row in changed_bins],
        "ct": [int(row.sum()) for row in changed_bins],
        "bp": [int(pos[0]), int(pos[-1])]}
    print(f"| {nm} | {dens[0].mean():.4f} | {dens[0].min():.4f} | {dens[0].max():.4f} | "
          f"**{np.abs(dens - dens[0]).max():.4f}** | "
          f"**{int((changed_bins[-1] > 0).sum())} of {NB}** |")
(REPO / "logs" / f"A_mismatch_map_chr21_seed{SEED}.json").write_text(
    json.dumps(out, separators=(",", ":")))
print(f"\n**Provenance.** chr21 · HG00096 diploid, half `eta_tau` per haplotype · "
      f"seed {SEED}, every `tau` drawn from the SAME uniform stream · common axis = the "
      f"{int(in_chain.sum()):,} positions both models can emit · density is PER CELL over "
      f"{int(known.sum()):,} target-known cells · `chain_span` PROVISIONAL.")
