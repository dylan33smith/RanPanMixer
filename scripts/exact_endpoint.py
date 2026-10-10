#!/usr/bin/env python
"""B-EVL-exact-endpoint: the primary endpoint, exactly, at a genome-wide share.

⚠ WHY THIS GATES THE `tau` GRID. Decision A makes the genome-wide utility the
PRIMARY pre-registered endpoint and requires it be reported EXACTLY rather than
estimated from draws. Decision B then gives each chromosome a SHARE of the one
budget, so every figure on record -- all measured with `beta` summing to 1 over
chr21 alone -- describes a larger tilt than the mechanism will actually apply.
`docs/plan.md` forbids pre-registering a grid the endpoint provably does not vary
over, so the grid cannot be pinned until the endpoint is re-measured at the scale
that will be used.

TWO ROUTES TO THE SAME QUANTITY, AND THEY ARE NOT EQUALLY STABLE:

  (1) MARGINALS -- `state_marginals` then `expected_mismatches`. A LINEAR
      functional of the posterior marginals: no differencing, no cancellation,
      and cross-checked against the toy enumerator to 4.8e-15.
  (2) DERIVATIVES -- central finite differences of `log Z_p` in `eta`, which give
      `E_Q[u]` and `Var_Q(u)` in three forward passes and are much cheaper.

⚠ Route (2) is the one every recorded `Var_Q(u)` came from, and its second
difference is a cancellation: the numerator is of order `Var * h^2` while each
`log Z_p` carries absolute error of order 1e-13. Shrinking `beta` to a
chromosome's share shrinks `Var` by the square of the share while leaving that
error where it is, so route (2) must degrade and this script measures where. If it
degrades inside the range we intend to use, the pre-registered endpoint has to name
route (1).

    python scripts/exact_endpoint.py --shares 1.0 0.01755 --json out.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mechanism import calibration as C  # noqa: E402
from mechanism import ffbs as F  # noqa: E402
from mechanism import hmm as H  # noqa: E402
from mechanism import utility as U  # noqa: E402

TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)

# chr21's candidate genome-wide shares, from `logs/B_map_coverage.md`. The
# allocation RULE is not decided -- see the residual noted in the report -- so the
# default sweep spans the plausible range rather than assuming one.
SHARE_BY_GENETIC_LENGTH = 62.2 / 3543.8      # 1.755%
SHARE_BY_TOTAL_BP = 46.709983 / 2875.0       # 1.624%
SHARE_UNIFORM = 1.0 / 22.0                   # 4.545%


def load(chrom: str = "chr21", target: str = "HG00096"):
    d = REPO / "data" / "cohort" / chrom
    hap = np.load(d / "haplotypes.npy")
    ic = np.load(d / "in_chain.npy").astype(bool)
    sup = np.load(d / "support.npy")
    gp = np.load(d / "genetic_pos.npy")
    pos = np.load(d / "positions.npy")
    m = H.from_arrays(hap[:, ic], gp[ic], pos[ic], cohort_id=chrom)
    tgt = np.load(REPO / "data" / "targets" / target / chrom / "path.npy")[ic]
    return m, hap[:, ic], sup, ic, tgt, F.baseline_log_norm(m)


def exact_mismatch(m, hap, tgt_1d, beta, eta, baseline) -> float:
    """Route (1): exact E[differ] for ONE chain, from the marginals."""
    psi = U.tilt_potentials(m.allowed, hap, tgt_1d, beta, eta)
    fwd = F.forward(m, psi)
    gamma = F.state_marginals(m, fwd, psi)
    return F.expected_mismatches(m, hap, tgt_1d, gamma)


def derivative_route(m, hap, tgt_1d, beta, eta, baseline, h):
    """Route (2): E_Q[u] and Var_Q(u) by finite differences of log Z_p.

    ⚠ Central where it can be, ONE-SIDED at `eta` < `h`. `eta_tau` is
    non-negative by construction, so a central stencil at `tau` = 0 asks for a
    negative exponent, which `tilt_potentials` correctly refuses. That is a real
    limitation of this route, not a detail: at the floor it is a forward
    difference and so first-order in `h` rather than second.
    """
    lo = eta - h
    onesided = lo < 0.0
    grid = (eta, eta + h, eta + 2 * h) if onesided else (lo, eta, eta + h)
    lz = []
    for e in grid:
        psi = U.tilt_potentials(m.allowed, hap, tgt_1d, beta, e)
        # the bracket is this chain's own reach; pass it rather than eta
        lz.append(F.log_z_p(F.forward(m, psi), baseline, e,
                            budget=e * float(beta.sum())))
    if onesided:
        mean = (lz[1] - lz[0]) / h
        var = (lz[2] - 2 * lz[1] + lz[0]) / h ** 2
    else:
        mean = (lz[2] - lz[0]) / (2 * h)
        var = (lz[2] - 2 * lz[1] + lz[0]) / h ** 2
    return mean, var, lz


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shares", type=float, nargs="+",
                    default=[1.0, SHARE_UNIFORM, SHARE_BY_GENETIC_LENGTH])
    ap.add_argument("--h", type=float, default=0.05)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    m, hap, sup, ic, tgt, baseline = load()
    known = (tgt[:, 0] >= 0) | (tgt[:, 1] >= 0)
    denom = int((tgt >= 0).sum())
    print(f"chr21 · T = {m.T:,} · K = {m.K} · HG00096 diploid · "
          f"{denom:,} target-known cells · baseline_log_norm {baseline:.6f}\n")

    out = {}
    for share in a.shares:
        util = U.beta_v(sup, in_chain=ic, share=share, chrom="chr21")
        label = f"share={share:.6g}"
        print(f"### {label}  (sum beta = {util.beta.sum():.6g})\n")
        print("| `tau` | `eta_tau` | exact E[differ] | shift vs `tau`=0 | "
              "share of known cells | E_Q[u] (deriv) | Var_Q(u) (deriv) |")
        print("|---|---|---|---|---|---|---|")
        rows, base_val = [], None
        for tau in TAUS:
            eta = C.eta_tau(tau)
            # diploid: each haplotype at eta/2, and E[differ] is the sum
            e_dif = sum(exact_mismatch(m, hap, tgt[:, hh], util.beta,
                                       eta / 2.0, baseline) for hh in (0, 1))
            if base_val is None:
                base_val = e_dif
            mean, var, lz = derivative_route(m, hap, tgt[:, 0], util.beta,
                                             eta / 2.0, baseline, a.h)
            rows.append({"tau": tau, "eta_tau": eta, "e_differ": e_dif,
                         "shift": e_dif - base_val, "e_q_u": mean,
                         "var_q_u": var, "log_z": lz[1]})
            print(f"| {tau} | {eta:.6f} | {e_dif:.6f} | {e_dif-base_val:+.6e} | "
                  f"{e_dif/denom:.6%} | {mean:.6e} | {var:.6e} |")
        out[label] = rows

        shifts = [r["shift"] for r in rows]
        mono = all(shifts[i + 1] <= shifts[i] + 1e-12 for i in range(len(shifts) - 1))
        print(f"\n- total movement across the grid: **{shifts[-1]:+.6e} cells**")
        print(f"- monotone non-increasing in `tau`: **{mono}**")
        # float64 headroom on the ENDPOINT, which is what makes route (1) usable
        rel = abs(shifts[-1]) / base_val if base_val else float("nan")
        print(f"- relative movement {rel:.3e} against float64 eps 2.22e-16 "
              f"-> **{rel/2.220446049250313e-16:.3e}x** the noise floor")
        vars_ = np.array([r["var_q_u"] for r in rows])
        print(f"- Var_Q(u) by derivative: min {vars_.min():.3e}, "
              f"max {vars_.max():.3e}, any non-positive: "
              f"**{bool((vars_ <= 0).any())}**\n")

    if a.json:
        a.json.write_text(json.dumps(out, indent=1))
        print(f"JSON -> {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
