#!/usr/bin/env python
"""Genome-wide `e_q_u` and `e_q_u_over_ceiling`, EXACTLY, with zero Monte-Carlo error.

⚠⚠ WHY THIS EXISTS. Every utility figure on record is chr21 at SHARE 1 -- the whole
privacy budget spent on one chromosome -- and reads `utility_retained` = 0.6931 at
`tau` = 0 through 0.50 and 0.6932 at `tau` = 0.9, i.e. the mechanism moves the fourth
decimal. At genome scope each chromosome takes its share (decision J), `mean beta_v`
drops from 3.27e-06 to 4.52e-08. Dylan asked for the measured number rather than an
extrapolation, and MEASURED 2026-10-12 (n = 1 target, HG00096, 22 autosomes, 44 exact
passes): `e_q_u` goes 7.089676460e-01 -> 7.089677294e-01 across `tau` = 0 to 0.99, a
movement of +8.345e-08, and `e_q_u_over_ceiling` is 0.858617 at BOTH ends.

⚠ TWO FIRST-ORDER PREDICTIONS OF THAT MOVEMENT BOTH FAILED, in OPPOSITE directions.
The measured shrink is a factor of **43.6** -- the current-build chr21-at-share-1
movement of 3.6368e-06 over the genome-wide 8.345e-08. ⚠ This read "against chr21's
2.94e-06 ... a factor of 35.2" until 2026-10-13; 2.94e-06 is the PRE-REBUILD figure from
the finite-difference route and dividing it by a current genome number is a cross-build
ratio. Scaling by `mean
beta_v` predicts 72; scaling by sum(beta_v^2) -- the correct first-order factor if
matches were independent, measured 9.290392e-06 on chr21 against 7.223240e-08
genome-wide -- predicts 128.6. So the movement is NOT a function of budget dilution
alone. The untested explanation is that allele matches are POSITIVELY CORRELATED along
the HMM copying path, so Var_Q(u) exceeds sum(beta_v^2) by a chromosome-dependent
factor; settling it needs Var_Q(u) from PAIR marginals, which this script does not
compute. Do not quote a mechanism for the 43.6 until that is measured.

⚠ THE FIGURE THAT MATTERS IS THE MOVEMENT, NOT THE LEVEL. `e_q_u_over_ceiling` is
0.8586 at `tau` = 0 too, where the release carries ZERO target information, because any
two humans share most of their alleles. Quoting the level as the mechanism's utility is the
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


def alt_totals(target: str) -> np.ndarray:
    """Per-haplotype count of chain positions carrying a NON-REF allele, over ALL 22
    autosomes -- the genome-wide normaliser for `--weights target_alt`.

    ⚠ Always all 22, whatever `--chroms` says: standing constraint 4 puts ONE budget on
    the genome, so a chromosome's non-ref positions take their share of the genome's,
    never 1 of their own. A `--chroms` subset just measures part of the genome-wide u.
    """
    tot = np.zeros(2, dtype=np.int64)
    for c in AUTOSOMES:
        ic = np.load(REPO / "data" / "cohort" / c / "in_chain.npy").astype(bool)
        p = np.load(REPO / "data" / "targets" / target / c / "path.npy")[ic]
        tot += (p > 0).sum(axis=0)
    assert (tot > 0).all(), f"{target}: a haplotype with no non-ref allele genome-wide"
    return tot


def rarity_weights(hap, tgt):
    """`w_v = -log f_v` for ONE haplotype on ONE chromosome's chain, else 0.

    `f_v` is the COHORT FREQUENCY of the target's own allele at `v`: carriers / called,
    over the donor panel. So an allele every donor carries gets weight -log(1) = 0 and
    contributes nothing, while a singleton gets -log(1/88) = 4.477. The weight is
    therefore BOUNDED by log(K), which matters -- an unbounded weight would break the
    [0,1] utility range that standing constraint 2 requires.

    ⚠ TWO CLASSES GET ZERO, both deliberately:
      * the target has no call -- nothing to match;
      * NO donor carries the target's allele (f = 0, so -log f is infinite) -- the
        release CANNOT match there at any eta, so weight spent on it is wasted by
        construction. `target_alt` does NOT do this and spends budget on unmatchable
        positions; that is one of the two ways this weighting differs from it.
    """
    carriers = (hap == tgt[None, :]).sum(axis=0)
    called = (hap >= 0).sum(axis=0)
    usable = (tgt >= 0) & (carriers > 0)
    f = np.ones_like(called, dtype=np.float64)
    np.divide(carriers, np.maximum(called, 1), out=f, where=usable)
    w = np.zeros(tgt.shape[0], dtype=np.float64)
    w[usable] = -np.log(f[usable])
    return w


def rarity_totals(target: str) -> np.ndarray:
    """Per-haplotype sum of `-log f_v` over ALL 22 autosomes, the genome-wide
    normaliser for `--weights target_rarity`. Always all 22, for the same reason
    `alt_totals` is: one budget on the genome (standing constraint 4).
    """
    tot = np.zeros(2, dtype=np.float64)
    for c in AUTOSOMES:
        d = REPO / "data" / "cohort" / c
        ic = np.load(d / "in_chain.npy").astype(bool)
        hap = np.load(d / "haplotypes.npy")[:, ic]
        p = np.load(REPO / "data" / "targets" / target / c / "path.npy")[ic]
        for k in (0, 1):
            tot[k] += float(rarity_weights(hap, p[:, k]).sum())
    assert (tot > 0).all(), f"{target}: a haplotype with zero rarity weight genome-wide"
    return tot


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
    ap.add_argument("--from-json", type=Path, default=None,
                    help="re-report a saved run instead of recomputing it")
    a = ap.parse_args(argv)

    if a.from_json:
        d = json.loads(a.from_json.read_text())
        rows = d["results"]
        taus = sorted({r["tau"] for r in rows})
        print(f"re-reporting {len(rows)} saved row(s) from {a.from_json}. "
              f"NO marginals computed.")
        _report(rows, taus, d["scope"])
        return 0

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
                        "e_q_u_over_ceiling": eu / ceil if ceil else float("nan"),
                        "chroms": len(a.chroms)})
            print(f"  {t} tau={tau:g}: E_Q[u] = {eu:.9e}, u_ceiling = {ceil:.9e}, "
                  f"e_q_u_over_ceiling = {eu / ceil:.6f}", flush=True)

    _report(out, a.taus, scope)
    if a.json:
        a.json.write_text(json.dumps({"scope": scope, "results": out}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


def _report(out, taus, scope) -> None:
    """Emit the table from `out`. Shared by a live run and `--from-json`, so a
    relabelled row never requires recomputing 44 exact passes."""
    n = len({r["target"] for r in out})
    print(f"\n{scope} exact utility. n = {n} target(s), ZERO Monte-Carlo error.\n")
    print("| metric | " + " | ".join(f"{x:g}" for x in taus) + " |")
    print("|" + "---|" * (len(taus) + 1))
    for name, key in (("`e_q_u` *", "e_q_u"), ("`u_ceiling`", "u_ceiling"),
                      ("`e_q_u_over_ceiling`", "e_q_u_over_ceiling")):
        cells = []
        for x in taus:
            vals = [r[key] for r in out if r["tau"] == x]
            assert vals, f"no saved rows at tau = {x}; the report would be silent"
            cells.append(f"{np.mean(vals):.9f}")
        print(f"| {name} | " + " | ".join(cells) + " |")
    if len(taus) > 1:
        lo = np.mean([r["e_q_u"] for r in out if r["tau"] == min(taus)])
        hi = np.mean([r["e_q_u"] for r in out if r["tau"] == max(taus)])
        print(f"\n⚠⚠ **THE MOVEMENT**: E_Q[u] goes {lo:.9e} -> {hi:.9e} across "
              f"tau = {min(taus):g} to {max(taus):g}, a change of **{hi - lo:+.3e}** "
              f"({100 * (hi - lo) / lo:+.6f}%).")
        # ⚠ the CURRENT-build chr21-at-share-1 movement, from
        # logs/B_beta_concentration.md:inv_support_all. The pre-rebuild 2.94e-06 that
        # stood here until 2026-10-13 made this a cross-build ratio.
        CHR21_SHARE1 = 3.6368e-06
        print(f"⚠ For comparison, chr21 at SHARE 1 moves {CHR21_SHARE1:.4e} (a factor of "
              f"{CHR21_SHARE1 / (hi - lo):.1f} more), and the toy at mean beta_v = 0.111 "
              f"moves 0.4502 -> 0.6477. The mechanism is not broken; the budget is "
              f"diluted.")
        print(f"⚠ That factor is NOT predicted by budget dilution: mean beta_v predicts "
              f"72, sum(beta_v^2) predicts 128.6. Positive correlation of matches along "
              f"the copying path is the UNTESTED explanation; it needs Var_Q(u) from pair "
              f"marginals, which this script does not compute.")


if __name__ == "__main__":
    sys.exit(main())
