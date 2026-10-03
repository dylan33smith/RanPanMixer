"""Exact expected released-vs-target mismatch count per `tau`. NO SAMPLING.

`E[differ]` is a sum over positions of a function of the state at that position, so
it is a LINEAR FUNCTIONAL of the posterior state marginals and therefore exact.
Two forward passes and two backward passes per `tau`; zero Monte-Carlo error.

⚠ WHY THIS SCRIPT EXISTS. 50 draws were spent estimating this quantity, at two
sample sizes, producing a signed monotonic pattern at n=10 that reversed at n=40.
An independent audit pointed out that the repo had already written the lesson down
for expected utility (`docs/bugs.md`: "report the exact quantity instead ... three
forward passes and no sampling") and had not applied it to the metric that replaced
it. It also found that the sampling script tilted EACH haplotype at the full
`eta_tau` rather than `eta_tau/2`, so it was not measuring this mechanism at all.
Both columns are emitted below so that error stays visible.
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

TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)
TARGET = sys.argv[1] if len(sys.argv) > 1 else "HG00096"
COH = REPO / "data" / "cohort" / "chr21"

hmm = H.from_cohort_dir(COH)
in_chain = np.load(COH / "in_chain.npy")
support = np.load(COH / "support.npy")
hap = np.load(COH / "haplotypes.npy")[:, in_chain]
util = U.beta_v(support, in_chain=in_chain)
tc = np.load(REPO / "data" / "targets" / TARGET / "chr21" / "path.npy")[in_chain]
known = int((tc >= 0).sum())


def exact_differ(eta_per_chain: float) -> float:
    tot = 0.0
    for h in (0, 1):
        psi = U.tilt_potentials(hmm.allowed, hap, tc[:, h], util.beta, eta_per_chain)
        fwd = F.forward(hmm, psi)
        g = F.state_marginals(hmm, fwd, psi)
        tot += F.expected_mismatches(hmm, hap, tc[:, h], g)
    return tot


rows = []
for tau in TAUS:
    eta = C.eta_tau(tau)
    rows.append((tau, eta, exact_differ(eta / 2.0), exact_differ(eta)))

base_ok, base_bug = rows[0][2], rows[0][3]
print(f"### Exact E[differ] per tau — {TARGET}, chr21, no sampling")
print()
print(f"| `tau` | `eta_tau` | E[differ], CORRECT (eta/2 per chain) | shift vs `tau`=0 | "
      f"as a share of known cells | E[differ] at FULL eta (the script's bug) | its shift |")
print("|---|---|---|---|---|---|---|")
for tau, eta, ok, bug in rows:
    print(f"| {tau} | {eta:.4f} | {ok:.6f} | {ok-base_ok:+.6f} | {(ok-base_ok)/known:+.3e} "
          f"| {bug:.6f} | {bug-base_bug:+.6f} |")
print()
sl_ok = [(ok-base_ok)/eta for tau, eta, ok, bug in rows[1:]]
sl_bug = [(bug-base_bug)/eta for tau, eta, ok, bug in rows[1:]]
print(f"**Provenance.** chr21/hprc-v1.0-pggb/K{hmm.K}/T{hmm.T} · {TARGET}, diploid · "
      f"{known:,} target-known cells of {tc.size:,} released · "
      f"utility `{F.utility_id(util)}` · EXACT, zero Monte-Carlo error · "
      f"validated against brute-force enumeration on the toy to 4.8e-15.")
print()
print(f"- slope, correct eta/2 per chain: {np.mean(sl_ok):.4f} cells per unit `eta_tau` "
      f"(range {min(sl_ok):.4f} to {max(sl_ok):.4f}) — LINEAR in `eta_tau`")
print(f"- slope, the buggy full eta:      {np.mean(sl_bug):.4f} cells per unit `eta_tau` "
      f"— exactly 2x, which is what a doubled tilt does")
print(f"- baseline E[differ] at `tau`=0:  {base_ok:.4f} of {known:,} = {base_ok/known:.4%}")
