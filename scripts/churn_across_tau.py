"""Does raising `tau` move a few positions, or reshuffle many and net to a few?

The exact `E[differ]` figures are NET. This measures the mismatch SET at a fixed seed
(common random numbers, the right paired design) and decomposes the change:

    fixed   -- mismatched at tau=0, matches now      (the tilt helping)
    broken  -- matched at tau=0, mismatches now      (the tilt hurting)
    net     -- broken - fixed, which is what E[differ] reports
    churn   -- fixed + broken, how many positions MOVED at all

⚠ BOTH MODELS RUN ON THE SAME POSITION AXIS. Model B cannot emit 1,222 positions (4
block units no donor covers), so both are restricted to the 304,664 both can emit.
Otherwise model B's mismatch count would be divided by the full chain's denominator,
which is the error the first version of this comparison made.
"""
from __future__ import annotations

import json
import os
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
SEEDS = (1000, 2000, 3000)   # ⚠ ONE target (HG00096). Seeds, not targets.
COH = REPO / "data" / "cohort" / "chr21"
# Outputs. Overridable with RANPANMIXER_SCRATCH (same convention as
# RANPANMIXER_VG); the default is repo-relative and gitignored, so a fresh
# checkout reruns without editing this file.
SCRATCH = Path(os.environ.get("RANPANMIXER_SCRATCH",
                              REPO / "logs" / "churn_across_tau"))
SCRATCH.mkdir(parents=True, exist_ok=True)

hap_full = np.load(COH / "haplotypes.npy")
gpos = np.load(COH / "genetic_pos.npy")
pos = np.load(COH / "positions.npy")
support = np.load(COH / "support.npy")
in_chain0 = np.load(COH / "in_chain.npy")
blocks = np.load(COH / "blocks_idx.npy")

# THE COMMON AXIS: drop what model B cannot emit, from BOTH models
dead = H.dead_block_positions(hap_full, blocks, in_chain0)
in_chain = in_chain0 & ~dead
hap = hap_full[:, in_chain]
util = U.beta_v(support, in_chain=in_chain)
tc = np.load(REPO / "data" / "targets" / "HG00096" / "chr21" / "path.npy")[in_chain]
known = tc >= 0
nk = int(known.sum())

models = {"C": H.from_arrays(hap_full, gpos, pos, in_chain=in_chain, blocks=None),
          "B": H.from_arrays(hap_full, gpos, pos, in_chain=in_chain, blocks=blocks)}

exact = {}
for nm, m in models.items():
    exact[nm] = {}
    for tau in TAUS:
        tot = 0.0
        for h in (0, 1):
            psi = U.tilt_potentials(m.allowed, hap, tc[:, h], util.beta,
                                    C.eta_tau(tau) / 2.0, unit_spans=m.unit_spans)
            tot += F.expected_mismatches(m, hap, tc[:, h],
                                         F.state_marginals(m, F.forward(m, psi), psi))
        exact[nm][tau] = tot

rows = []
for nm, m in models.items():
    for seed in SEEDS:
        mism = {}
        for tau in TAUS:
            r = F.release(m, hap, tc, util, tau, seed=seed, target_id="HG00096")
            mism[tau] = (r.released_path != tc) & known
        np.save(SCRATCH / f"mism_{nm}_{seed}.npy",
                np.stack([mism[t] for t in TAUS]))
        base = mism[0.0]
        prev = base
        for tau in TAUS:
            cur = mism[tau]
            rows.append(dict(
                model=nm, seed=seed, tau=tau, n=int(cur.sum()),
                fixed0=int((base & ~cur).sum()), broken0=int((~base & cur).sum()),
                fixedp=int((prev & ~cur).sum()), brokenp=int((~prev & cur).sum())))
            prev = cur

json.dump({"rows": rows, "taus": list(TAUS), "seeds": list(SEEDS),
           "exact": {k: {str(t): v for t, v in d.items()} for k, d in exact.items()},
           "known_cells": nk, "positions": int(in_chain.sum())},
          open(SCRATCH / "churn.json", "w"), indent=1)

print(f"### Mismatch-set churn across `tau`, same seed — common axis, "
      f"{int(in_chain.sum()):,} positions / {nk:,} known cells\n")
for nm in models:
    print(f"#### Model {nm} — every cell is a MEAN over {len(SEEDS)} seeds, so fractions are means of integers")
    print(f"| `tau` | mismatches (mean of {len(SEEDS)}) | net vs `tau`=0 | fixed | broken "
          f"| **CHURN (mean)** | churn per seed | churn vs previous |")
    print("|---|---|---|---|---|---|---|---|")
    for tau in TAUS:
        rs = sorted([r for r in rows if r["model"] == nm and r["tau"] == tau],
                    key=lambda r: r["seed"])
        n = np.mean([r["n"] for r in rs]); fx = np.mean([r["fixed0"] for r in rs])
        bk = np.mean([r["broken0"] for r in rs]); cp = np.mean([r["fixedp"] + r["brokenp"] for r in rs])
        per = " / ".join(str(r["fixed0"] + r["broken0"]) for r in rs)
        print(f"| {tau} | {n:,.1f} | {bk-fx:+,.1f} | {fx:,.1f} | {bk:,.1f} | "
              f"**{fx+bk:,.1f}** | {per} | {cp:,.1f} |")
    print(f"\nexact E[differ]: " + " · ".join(f"`tau`={t}: {exact[nm][t]:,.3f}" for t in TAUS))
    print()
print(f"**Provenance.** chr21 · HG00096 diploid, half `eta_tau` per haplotype · "
      f"ONE target, {len(SEEDS)} seeds {SEEDS} (NOT 10 targets) · each column within a "
      f"seed uses the SAME "
      f"uniform stream (paired) · common axis = the 304,664 positions BOTH models can "
      f"emit · utility `{F.utility_id(util)}` · `chain_span` PROVISIONAL.")
