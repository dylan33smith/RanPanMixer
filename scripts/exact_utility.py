#!/usr/bin/env python
"""Genome-wide E_Q[u] and `utility_retained`, EXACTLY, with zero Monte-Carlo error.

⚠⚠ WHY THIS EXISTS. Every utility figure on record is chr21 at SHARE 1 -- the whole
privacy budget spent on one chromosome -- and reads `utility_retained` = 0.6931 at
`tau` = 0 through 0.50 and 0.6932 at `tau` = 0.9, i.e. the mechanism moves the fourth
decimal. At genome scope each chromosome takes its share (decision J), `mean beta_v`
drops from 3.27e-06 to 4.52e-08, and the movement should shrink by about the same
factor of 72. Dylan asked for the measured number rather than that extrapolation.

⚠ THE FIGURE THAT MATTERS IS THE MOVEMENT, NOT THE LEVEL. `utility_retained` is ~0.69
at `tau` = 0 too, where the release carries ZERO target information, because any two
humans share most of their alleles. Quoting the level as the mechanism's utility is the
free-lunch error `target_fidelity` was defined to prevent. So this reports
E_Q[u] - E_R[u] as its own row.

ROUTE. `E_Q[u]` is a LINEAR functional of the posterior state marginals, so one
forward + one `state_marginals` pass per (chromosome, haplotype, tau) gives it with no
sampling and no differencing. `tau` = 0 needs the UNTILTED marginals, which is the same
pass with all-ones potentials.

    python scripts/exact_utility.py --targets 1 --taus 0 0.99
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

from mechanism import calibration as C                    # noqa: E402
from mechanism import ffbs, hmm, utility                  # noqa: E402

AUTOSOMES = [f"chr{c}" for c in range(1, 23)]


def shares() -> dict:
    """`beta_share` by variant count, `T_c / T_total` (decision J)."""
    T = {c: int(np.load(REPO / "data" / "cohort" / c / "in_chain.npy").sum())
         for c in AUTOSOMES}
    tot = sum(T.values())
    return {c: T[c] / tot for c in AUTOSOMES}, T, tot


def chrom_eu(chrom: str, target: str, tau: float, share: float) -> tuple:
    """`E_Q[u]` and `u_ceiling` for one chromosome at its share. Exact."""
    d = REPO / "data" / "cohort" / chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=chrom)
    util = utility.beta_v(np.load(d / "support.npy"), in_chain=ic,
                          share=share, chrom=chrom)
    p = np.load(REPO / "data" / "targets" / target / chrom / "path.npy")[ic]
    eta = C.eta_tau(tau)
    base = ffbs.baseline_log_norm(m)
    e_u = 0.0
    for col in (0, 1):
        tgt = p[:, col]
        # ⚠ at eta = 0 the potentials are all-ones and this is the BASELINE posterior,
        # which is exactly what the tau = 0 floor means -- no special case needed.
        psi = utility.tilt_potentials(m.allowed, hap, tgt, util.beta, eta / 2.0)
        fwd = ffbs.forward(m, psi)
        gam = ffbs.state_marginals(m, fwd, psi)
        known = tgt >= 0
        match = (hap == tgt[None, :]) & known[None, :]
        # ⚠ 0.5 * each haplotype: `diploid_utility` is the MEAN, matching the density
        # exp(eta * u_mean) that the sampler draws from.
        e_u += 0.5 * float((gam * match.T).sum(axis=1) @ util.beta)
    ceil = utility.diploid_utility(utility.u_ceiling(p[:, 0], util.beta),
                                   utility.u_ceiling(p[:, 1], util.beta))
    return e_u, ceil


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=int, default=1)
    ap.add_argument("--taus", type=float, nargs="+", default=[0.0, 0.99])
    ap.add_argument("--chroms", nargs="+", default=AUTOSOMES)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    sh, T, tot = shares()
    sub = {c: sh[c] for c in a.chroms}
    scope = "GENOME-WIDE" if len(a.chroms) == 22 else f"{len(a.chroms)} chromosome(s)"
    print(f"{scope}: T = {sum(T[c] for c in a.chroms):,} of {tot:,}, "
          f"share covered = {sum(sub.values()):.6f}")
    print(f"⚠ The MOVEMENT is the figure, not the level: `utility_retained` is ~0.69 at "
          f"tau = 0 too, where the release carries ZERO target information.\n")

    tids = sorted(p.name for p in (REPO / "data" / "targets").iterdir()
                  if p.is_dir())[: a.targets]
    out = []
    t0 = time.time()
    for t in tids:
        for tau in a.taus:
            eu = ceil = 0.0
            for c in a.chroms:
                e, k = chrom_eu(c, t, tau, sh[c])
                eu += e
                ceil += k
                print(f"    {t} tau={tau:g} {c}: E_Q[u]={e:.6e} "
                      f"({time.time() - t0:.0f}s)", flush=True)
            out.append({"target": t, "tau": tau, "e_q_u": eu, "u_ceiling": ceil,
                        "utility_retained": eu / ceil if ceil else float("nan"),
                        "chroms": len(a.chroms)})
            print(f"  {t} tau={tau:g}: E_Q[u] = {eu:.9e}, u_ceiling = {ceil:.9e}, "
                  f"utility_retained = {eu / ceil:.6f}", flush=True)

    print(f"\n{scope} exact utility. n = {len(tids)} target(s), ZERO Monte-Carlo error.\n")
    print("| metric | " + " | ".join(f"{x:g}" for x in a.taus) + " |")
    print("|" + "---|" * (len(a.taus) + 1))
    for name, key in (("`e_q_u`", "e_q_u"), ("`u_ceiling`", "u_ceiling"),
                      ("`utility_retained` *", "utility_retained")):
        cells = [f"{np.mean([r[key] for r in out if r['tau'] == x]):.9f}"
                 for x in a.taus]
        print(f"| {name} | " + " | ".join(cells) + " |")
    if len(a.taus) > 1:
        lo = np.mean([r["e_q_u"] for r in out if r["tau"] == min(a.taus)])
        hi = np.mean([r["e_q_u"] for r in out if r["tau"] == max(a.taus)])
        print(f"\n⚠⚠ **THE MOVEMENT**: E_Q[u] goes {lo:.9e} -> {hi:.9e} across "
              f"tau = {min(a.taus):g} to {max(a.taus):g}, a change of **{hi - lo:+.3e}** "
              f"({100 * (hi - lo) / lo:+.6f}%).")
        print(f"⚠ For comparison, chr21 at SHARE 1 moves 2.94e-06, and the toy at "
              f"mean beta_v = 0.111 moves 0.4502 -> 0.6477. The mechanism is not broken; "
              f"the budget is diluted.")
    if a.json:
        a.json.write_text(json.dumps({"scope": scope, "results": out}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
