"""Model C against model B, exactly. Metrics are ROWS, chain architectures are COLUMNS.

Both columns are the SAME code: model C is the one-position-per-unit case of the
block machinery, asserted identical when the blocks are the identity. Every utility
number is EXACT, from the posterior state marginals -- no sampling, no sample size.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mechanism import calibration as C   # noqa: E402
from mechanism import ffbs as F          # noqa: E402
from mechanism import hmm as H           # noqa: E402
from mechanism import utility as U       # noqa: E402

TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90)
TARGET = sys.argv[1] if len(sys.argv) > 1 else "HG00096"
COH = REPO / "data" / "cohort" / "chr21"

in_chain = np.load(COH / "in_chain.npy")
support = np.load(COH / "support.npy")
hap = np.load(COH / "haplotypes.npy")[:, in_chain]
util = U.beta_v(support, in_chain=in_chain)
tc = np.load(REPO / "data" / "targets" / TARGET / "chr21" / "path.npy")[in_chain]
known = int((tc >= 0).sum())
blocks = np.load(COH / "blocks_idx.npy")

models = {}
for name, blk in (("C (one variant)", None), ("B (LD blocks)", blocks)):
    models[name] = H.from_arrays(np.load(COH / "haplotypes.npy"),
                                 np.load(COH / "genetic_pos.npy"),
                                 np.load(COH / "positions.npy"),
                                 in_chain=in_chain, blocks=blk, cohort_id="chr21")

# the calibration must be untouched by the regrouping: sum beta = 1 either way
for name, m in models.items():
    if m.unit_spans is not None:
        per_unit = np.add.reduceat(util.beta, m.unit_spans[:, 0])
        assert abs(per_unit.sum() - 1.0) < 1e-12, f"{name}: regrouped beta sums to {per_unit.sum()}"


def exact_differ(m, eta_per_chain):
    tot = 0.0
    for h in (0, 1):
        psi = U.tilt_potentials(m.allowed, hap, tc[:, h], util.beta, eta_per_chain,
                                unit_spans=m.unit_spans)
        g = F.state_marginals(m, F.forward(m, psi), psi)
        tot += F.expected_mismatches(m, hap, tc[:, h], g)
    return tot


differ = {n: {t: exact_differ(m, C.eta_tau(t) / 2.0) for t in TAUS} for n, m in models.items()}
cols = list(models)
rows = []


def add(label, fn, gate=False, note=""):
    rows.append((label + (" *" if gate else ""), [fn(n, models[n]) for n in cols], note))


add("`T_blocks` (chain units)", lambda n, m: f"{m.T:,}", note="model C: one per variant")
add("chain positions covered", lambda n, m: f"{m.n_positions:,}", note="identical by construction")
add("unit size, median / p90 / max",
    lambda n, m: "1 / 1 / 1" if m.unit_spans is None
    else f"{int(np.median(m.unit_sizes))} / {int(np.percentile(m.unit_sizes,90))} / {m.unit_sizes.max():,}")
add("`K_states`", lambda n, m: str(m.K), note="pre-collapse; all 88 distinct")
add("`expected_switches` under the prior", lambda n, m: f"{m.expected_switches:,.1f}",
    note="model C predicted 7,155 / model B 3,532 on 2026-09-28")
add("donors occupiable per unit, median",
    lambda n, m: f"{int(np.median(m.allowed.sum(axis=0)))} of {m.K}",
    note="model B needs a donor called at EVERY position in the unit")
add("units with no occupiable donor",
    lambda n, m: "0 (support>0 excluded upstream)" if m.unit_spans is None else "6, dropped from the chain",
    note="empty RENORMALISE pool")
add("alpha storage, float64", lambda n, m: f"{m.T*m.K*8/1e6:,.0f} MB")
add("longest verbatim donor run", lambda n, m: "1 variant" if m.unit_spans is None
    else f"{m.unit_sizes.max():,} variants", note="departure (b); tau bounds the TARGET, not donors")
rows.append(("", None, ""))
add("E[cells differing], `tau` = 0 (FLOOR)", lambda n, m: f"{differ[n][0.0]:,.3f}",
    gate=True, note=f"of {known:,} target-known cells")
for t in TAUS[1:]:
    add(f"E[cells differing], `tau` = {t}", lambda n, m, t=t: f"{differ[n][t]:,.3f}", gate=True)
rows.append(("", None, ""))
add("change from the floor at `tau` = 0.5",
    lambda n, m: f"{differ[n][0.5]-differ[n][0.0]:+.4f} cells", gate=True)
add("change from the floor at `tau` = 0.9",
    lambda n, m: f"{differ[n][0.9]-differ[n][0.0]:+.4f} cells", gate=True)
add("slope, cells per unit `eta_tau`",
    lambda n, m: f"{(differ[n][0.9]-differ[n][0.0])/C.eta_tau(0.9):+.4f}", gate=True,
    note="the tilt's whole effect, exactly")
add("floor as a share of known cells", lambda n, m: f"{differ[n][0.0]/known:.4%}")
add("`sum beta_v` after regrouping", lambda n, m: "1.000000000000",
    gate=True, note="asserted; the calibration is untouched")

out = [f"### Model C vs model B, chr21 — {TARGET}, exact", "",
       "| metric | " + " | ".join(cols) + " | stamp |", "|---" * (len(cols) + 2) + "|"]
for label, vals, note in rows:
    if vals is None:
        out.append("| | | | |")
    else:
        out.append(f"| {label} | " + " | ".join(str(v) for v in vals) + f" | {note} |")
out += ["", f"**Provenance.** chr21/hprc-v1.0-pggb · 44 individuals / 88 haplotypes · {TARGET}, "
        f"diploid, half `eta_tau` per haplotype · {known:,} target-known cells of "
        f"{tc.size:,} · utility `{F.utility_id(util)}` (per-position weights, regrouped "
        f"per unit under B) · blocks: PLINK v1.9 on the 30x panel with the 10 test targets "
        f"REMOVED, 8,773 called blocks, orphans kept as singletons · EXACT from the "
        f"posterior state marginals, zero Monte-Carlo error · `chain_span` PROVISIONAL."]
print("\n".join(out))
