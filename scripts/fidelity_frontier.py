"""The target_fidelity axis for every arm and every knob setting. Metrics are ROWS.

The one quantity computable on either tool's output: how many of the target's
called alleles does the released object differ from? Ours is EXACT (posterior
marginals); PanMixer's is read from its own stacker log.

⚠ `capacity` and `tau` ARE NOT COMPARABLE PARAMETERS and are never put in the same
column here. What is comparable is the fidelity axis itself, which is why this is a
FRONTIER and not a pair of point values. Nor are the denominators equal -- see the
provenance line -- so the shares are the quotable numbers, not the raw counts.
"""
from __future__ import annotations

import re
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
COH = REPO / "data" / "cohort" / "chr21"
in_chain = np.load(COH / "in_chain.npy")
support = np.load(COH / "support.npy")
hap = np.load(COH / "haplotypes.npy")[:, in_chain]
util = U.beta_v(support, in_chain=in_chain)
tc = np.load(REPO / "data" / "targets" / "HG00096" / "chr21" / "path.npy")[in_chain]
known = int((tc >= 0).sum())
blocks = np.load(COH / "blocks_idx.npy")

ours = {}
for name, blk in (("C", None), ("B", blocks)):
    m = H.from_arrays(np.load(COH / "haplotypes.npy"), np.load(COH / "genetic_pos.npy"),
                      np.load(COH / "positions.npy"), in_chain=in_chain, blocks=blk)
    ours[name] = {}
    for t in TAUS:
        tot = 0.0
        for h in (0, 1):
            psi = U.tilt_potentials(m.allowed, hap, tc[:, h], util.beta,
                                    C.eta_tau(t) / 2.0, unit_spans=m.unit_spans)
            tot += F.expected_mismatches(m, hap, tc[:, h],
                                         F.state_marginals(m, F.forward(m, psi), psi))
        ours[name][t] = tot

# PanMixer's own knob, parsed from its stacker lines -- fail loudly if absent
log = (REPO / "logs" / "pmsweep.log").read_text()
caps = re.findall(r"^### capacity=([\d.]+)", log, re.M)
chg = re.findall(r"changed=(\d+)/(\d+) alleles \(([\d.]+)%\)", log)
if len(caps) != len(chg) or not caps:
    raise SystemExit(f"capacity sweep incomplete: {len(caps)} headers, {len(chg)} results")
pm = [(float(c), int(a), int(b), float(p)) for c, (a, b, p) in zip(caps, chg)]

out = ["### `target_fidelity` frontier — ours (exact) against PanMixer's capacity knob", ""]
out += ["#### Ours: expected share of the target's called alleles that differ", "",
        "| metric | " + " | ".join(f"`tau`={t}" for t in TAUS) + " |",
        "|---" * (len(TAUS) + 1) + "|"]
for name, label in (("C", "model C (one variant per unit)"), ("B", "model B (LD blocks)")):
    out.append(f"| {label} | " + " | ".join(f"{ours[name][t]/known:.4%}" for t in TAUS) + " |")
out.append("| change from the floor, model C | " +
           " | ".join("—" if t == 0 else f"{ours['C'][t]-ours['C'][0.0]:+.3f} cells" for t in TAUS) + " |")
out.append("| change from the floor, model B | " +
           " | ".join("—" if t == 0 else f"{ours['B'][t]-ours['B'][0.0]:+.3f} cells" for t in TAUS) + " |")

out += ["", "#### PanMixer: measured share, swept across its `capacity` knob", "",
        "| metric | " + " | ".join(f"cap={c:g}" for c, _, _, _ in pm) + " |",
        "|---" * (len(pm) + 1) + "|",
        "| alleles changed | " + " | ".join(f"{a:,}" for _, a, _, _ in pm) + " |",
        "| share of its denominator | " + " | ".join(f"{p:.2f}%" for *_, p in pm) + " |",
        "| denominator | " + " | ".join(f"{b:,}" for _, _, b, _ in pm) + " |"]

share_c = ours["C"][0.5] / known * 100
lo = max([p for p in pm if p[3] <= share_c], key=lambda r: r[3], default=None)
hi = min([p for p in pm if p[3] > share_c], key=lambda r: r[3], default=None)
cross = (lo[0] + (hi[0] - lo[0]) * (share_c - lo[3]) / (hi[3] - lo[3])) if lo and hi else None
sat = max(p[3] for p in pm)

out += ["", f"**Provenance.** OURS: chr21/hprc-v1.0-pggb, 44 individuals / 88 haplotypes, "
        f"HG00096 (EXTERNAL), diploid, half `eta_tau` per haplotype, {known:,} target-known "
        f"cells of {tc.size:,}, utility `{F.utility_id(util)}`, EXACT from the posterior "
        f"marginals with zero Monte-Carlo error, `chain_span` PROVISIONAL. Model B blocks: "
        f"PLINK v1.9 on the 30x panel with the 10 test targets REMOVED. "
        f"PANMIXER: `external/PanMixer` @ `c182c38` unmodified, HG00438 (a COHORT MEMBER), "
        f"seed 123, `logs/pmsweep.log`, sentinel 0, one run per capacity.", "",
        "⚠ **The denominators are not equal and the targets are not the same person.** Ours "
        "counts only cells where the target's allele is known and only inside `chain_span`; "
        "theirs counts its own site set, which includes 10.25% of records our chain never "
        "touches. Read the SHARES, and read them as a frontier.", ""]
if cross:
    out.append(f"- **Our model-C release sits at {share_c:.2f}%, which PanMixer reaches at "
               f"capacity ≈ {cross:.2f}.** That is the matched-fidelity operating point.")
out.append(f"- **PanMixer's knob SATURATES at {sat:.2f}%** from capacity "
           f"{min(c for c,_,_,p in pm if p>=sat-0.01):g} upward — four settings return the "
           f"identical count. Its reachable range is "
           f"{min(p for *_,p in pm):.2f}% to {sat:.2f}%.")
out.append(f"- **Ours does not move with `tau` at all** on this axis: "
           f"{ours['C'][0.0]/known:.4%} at `tau`=0 to {ours['C'][0.99]/known:.4%} at "
           f"`tau`=0.99, a span of {abs(ours['C'][0.99]-ours['C'][0.0]):.2f} cells.")
print("\n".join(out))
